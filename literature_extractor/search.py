"""Search public article metadata and export a candidate reading list."""

from concurrent.futures import ThreadPoolExecutor
from html import unescape
from threading import Lock
from time import monotonic, sleep
from urllib.error import HTTPError
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen
import argparse
import csv
import json
import re
import sys

EUROPE = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search'
CROSSREF = 'https://api.crossref.org/works'
CATEGORIES = {'synthesis': '全合成', 'assignment': '专门归属', 'revision': '归属修订', 'isolation': '分离与结构鉴定'}
CROSSREF_GATE = Lock()
LAST_CROSSREF = 0.0


def plain(value):
    if not isinstance(value, str):
        return ''
    decoded = re.sub(r'</?(?:i|b|sub|sup)\b[^>]*>', '', unescape(value), flags=re.I)
    return ' '.join(re.sub(r'<[^>]*>', ' ', decoded).split())


def read(url, params):
    if url == CROSSREF:
        return read_crossref(params)
    return request_json(url, params)


def request_json(url, params):
    req = Request(url + '?' + urlencode(params), headers={'User-Agent': 'AtomNMR-Literature/1.0 (local metadata client)',
                                                         'Accept': 'application/json'})
    with urlopen(req, timeout=12) as response:
        return json.load(response)


def read_crossref(params):
    global LAST_CROSSREF
    with CROSSREF_GATE:
        sleep(max(0, 1.1 - (monotonic() - LAST_CROSSREF)))
        for attempt in range(2):
            LAST_CROSSREF = monotonic()
            try:
                return request_json(CROSSREF, params)
            except HTTPError as exc:
                if exc.code != 429 or attempt:
                    raise
                wait = exc.headers.get('Retry-After', '')
                if wait.isdigit() and int(wait) > 10:
                    raise
                sleep(max(2, int(wait) if wait.isdigit() else 2))


def queries(name):
    named = f'"{name}" AND ' if name else ''
    return {
        'synthesis': (f'{named}"total synthesis" AND "NMR assignments" AND "natural product"',
                      f'{name} total synthesis NMR assignments natural product', 'total synthesis'),
        'assignment': (f'{named}("complete NMR assignment" OR "chemical shift assignments" OR "NMR assignments")',
                       f'{name} complete NMR assignment chemical shift assignments', f'{name} NMR assignments'),
        'revision': (f'({named}"revised NMR assignments") OR ({named}"structure revision" AND "HMBC")',
                     f'{name} revised NMR assignments structure revision HMBC', f'{name} structure revision'),
        'isolation': (f'{named}"isolation" AND "structure elucidation" AND "NMR"',
                      f'{name} isolation structure elucidation NMR', f'{name} isolation structure elucidation'),
    }


def fetch(provider, category, terms, size):
    europe, broad, title = terms
    if provider == 'Europe PMC':
        raw = read(EUROPE, {'query': europe, 'resultType': 'core', 'format': 'json', 'pageSize': size})
        for row in raw.get('resultList', {}).get('result', []):
            pmcid = row.get('pmcid') or ''
            doi = (row.get('doi') or '').lower()
            link = doi_url(doi) if doi else f'https://europepmc.org/article/{row["source"]}/{row["id"]}'
            open_link = f'https://europepmc.org/articles/{pmcid}' if pmcid else ''
            yield dict(title=plain(row.get('title')), doi=doi, journal=plain(row.get('journalTitle') or
                (row.get('journalInfo') or {}).get('journal', {}).get('title')), year=str(row.get('pubYear') or ''),
                url=link, open_url=open_link, abstract=plain(row.get('abstractText')), provider=provider)
    else:
        params = {'query.bibliographic': broad, 'query.title': title, 'filter': 'type:journal-article', 'rows': size}
        for row in read(CROSSREF, params).get('message', {}).get('items', []):
            doi = row['DOI'].lower()
            date = (row.get('published') or {}).get('date-parts', [[]])[0]
            yield dict(title=plain((row.get('title') or [''])[0]), doi=doi, url=doi_url(doi), open_url='',
                       journal=plain((row.get('container-title') or [''])[0]), year=str(date[0]) if date else '',
                       abstract=plain(row.get('abstract')), provider=provider)


def doi_url(doi):
    return 'https://doi.org/' + quote(doi, safe='/;().:-_')


def relevance(row, category, name):
    title, abstract = row['title'].casefold(), row['abstract'].casefold()
    whole = title + ' ' + abstract
    if name and name.casefold() not in whole:
        return 0, ''
    if re.search(r'\b(?:backbone and side.chain|(?:kunitz|protein|binding|sh[23]) domain|disordered proteins?|'
                 r'protein resonance|ribosomal protein|rna|dna)\b', title):
        return 0, ''
    assignment = bool(re.search(r'(?:assign.{0,38}(?:nmr|chemical.shift|spectr|signal)|(?:nmr|chemical.shift).{0,40}assign)', whole))
    nmr = 'nmr' in whole or 'nuclear magnetic resonance' in whole
    two_d = bool(re.search(r'\b(?:hsqc|hmqc|hmbc|cosy)\b', whole))
    title_patterns = {'synthesis': 'total synthes',
                      'assignment': r'assign|nmr|magnetic resonance',
                      'revision': r'revis|reassign|structur.{0,28}(?:correct|elucid|confirm)',
                      'isolation': r'isolat|new.{0,40}(?:product|alkaloid|terpen|polyketide|peptide|metabolite)'}
    title_match = bool(re.search(title_patterns[category], title))
    match = title_match or (bool(abstract) and bool(re.search(title_patterns[category], abstract)))
    if not match or (category == 'assignment' and not assignment) or (not title_match and not nmr):
        return 0, ''
    info = '题名或摘要提到 NMR 归属' if assignment else '题名或摘要提到二维 NMR' if two_d else '候选论文，归属表需打开确认'
    return 1 + 2*title_match + 3*assignment + 2*two_d + nmr, info


def search(name='', limit=8, categories=tuple(CATEGORIES)):
    name = re.sub(r'["\\\r\n]', ' ', name).strip()
    terms = queries(name)
    selected = {category: [] for category in categories}
    warnings = set()
    with ThreadPoolExecutor(max_workers=5) as pool:
        jobs = {pool.submit(lambda p, c: list(fetch(p, c, terms[c], min(75, limit*6))), provider, category):
                (provider, category) for provider in ('Crossref', 'Europe PMC') for category in categories}
        for job in jobs:
            provider, category = jobs[job]
            try:
                rows = job.result()
            except (OSError, ValueError, KeyError) as exc:
                warnings.add(f'{provider} · {CATEGORIES[category]} · {exc}')
                continue
            for rank, row in enumerate(rows):
                score, evidence = relevance(row, category, name)
                if row['title'] and score:
                    selected[category].append((score, rank, {**row, 'evidence': evidence}))
    results, by_key, title_keys = [], {}, {}
    for category, choices in selected.items():
        count, seen = 0, set()
        for _, _, row in sorted(choices, key=lambda x: (-x[0], x[1], x[2]['doi'])):
            title_key = re.sub(r'\W+', '', row['title'].casefold())
            key = row['doi'] or title_keys.get(title_key) or title_key
            title_keys[title_key] = key
            if key in seen:
                existing = by_key[key]
                existing['open_url'] = existing['open_url'] or row['open_url']
                continue
            if count >= limit:
                continue
            seen.add(key)
            count += 1
            if key not in by_key:
                entry = {k: v for k, v in row.items() if k != 'abstract'}
                entry['categories'] = []
                by_key[key] = entry
                results.append(entry)
            by_key[key]['categories'].append(CATEGORIES[category])
            by_key[key]['open_url'] = by_key[key]['open_url'] or row['open_url']
    return {'results': results, 'warnings': sorted(warnings), 'queries': {CATEGORIES[k]: terms[k][0] for k in selected}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--name', default='')
    parser.add_argument('--limit', type=int, choices=range(1, 21), default=8)
    parser.add_argument('--category', action='append', choices=CATEGORIES)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    report = search(args.name, args.limit, args.category or tuple(CATEGORIES))
    fields = ('title', 'url', 'open_url', 'doi', 'journal', 'year', 'categories', 'evidence', 'provider')
    with open(args.output, 'w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for item in report['results']:
            row = {**item, 'categories': ' / '.join(item['categories'])}
            writer.writerow({k: ("'" + v if v[:1] in ('=', '+', '-', '@') else v) for k, v in row.items()})
    for warning in report['warnings']:
        print(warning, file=sys.stderr)
    print(f'{len(report["results"])} results: {args.output}')


if __name__ == '__main__':
    main()
