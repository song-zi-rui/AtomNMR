"""Read text-based NMR assignment tables while preserving literal evidence."""

import json
import re
import sys
import unicodedata
from pathlib import Path


NUMBER = r'[+−-]?\d+(?:\.\d+)?'
SHIFT = re.compile(rf'({NUMBER})(?:[–—−-]({NUMBER}))?')
DOI = re.compile(r'10\.\d{4,9}/[-._;()/:A-Z0-9]+', re.I)


def normalize(value):
    return unicodedata.normalize('NFKC', value).casefold().replace('δ', 'd')


def headings(words):
    for carbon in words:
        if not re.fullmatch(r'd\(?13c\)?(?:\[ppm\])?', normalize(carbon['text'])):
            continue
        row = sorted((w for w in words if abs(w['top'] - carbon['top']) < 2.8), key=lambda w: w['x0'])
        hydrogen = next((w for w in row if w['x0'] < carbon['x0']
                         and re.fullmatch(r'd\(?1h\)?(?:\[ppm\])?', normalize(w['text']))), None)
        if not hydrogen:
            continue
        left = [w for w in row if hydrogen['x0'] < w['x0'] < carbon['x0']]
        right = [w for w in row if carbon['x0'] < w['x0'] < carbon['x0'] + 130]
        definitions = [
            ('h', hydrogen), ('mult', next((w for w in left if normalize(w['text']) in ('multiplicity', 'mult.')), None)),
            ('nh', next((w for w in left if normalize(w['text']) in ('no.', 'nh', 'integral')), None)),
            ('ha', next((w for w in left if normalize(w['text']).startswith('assignment')), None)),
            ('j', next((w for w in left if re.fullmatch(r'j(?:\[hz\])?', normalize(w['text']))), None)),
            ('c', carbon), ('ca', next((w for w in right if normalize(w['text']).startswith('assignment')), None)),
        ]
        if definitions[3][1] and definitions[6][1]:
            yield [(k, w) for k, w in definitions if w]


def parse_peak(raw, assignment, nucleus, row, cells):
    match = SHIFT.fullmatch(raw)
    if not match:
        raise ValueError(f'无法识别位移 {raw}')
    values = [float(part.replace('−', '-')) for part in match.groups() if part is not None]
    letter = 'H' if nucleus == '1H' else 'C'
    marks = re.findall(r'\[([a-z])\]', assignment, re.I)
    labels = []
    source_label = re.sub(r'\[[a-z]\]', '', assignment, flags=re.I).replace('’', '′').replace("'", '′')
    source_label = source_label.strip(',')
    for token in re.split(r'[,;]', source_label) if source_label and source_label not in ('–', '-','—') else []:
        atom = re.fullmatch(rf'{letter}(\d+)([′a-z])?(?:[–-]{letter}?(\d+))?', token)
        if not atom:
            raise ValueError(f'无法识别原子编号 {token}')
        first, suffix, last = atom.groups()
        if last:
            if suffix or not 0 <= int(last) - int(first) < 100:
                raise ValueError(f'无法识别编号范围 {token}')
            labels.extend(f'{letter}{n}' for n in range(int(first), int(last) + 1))
        else:
            labels.append(f'{letter}{first}{suffix or ""}')
    nh = cells.get('nh', '').removesuffix('H') if nucleus == '1H' else ''
    j = cells.get('j', '') if nucleus == '1H' else ''
    return {
        'row': row, 'shift_text': raw, 'shift_ppm': values[0] if len(values) == 1 else None,
        'range_ppm': sorted(values) if len(values) > 1 else None,
        'source_assignment_text': assignment, 'source_labels': labels, 'source_footnotes': marks,
        'multiplicity': (cells.get('mult') or None) if nucleus == '1H' else None,
        'reported_nH': float(nh) if nh and nh not in ('–', '-', '—') else None,
        'j_values_hz': [float(part) for part in j.split(',')] if j and j not in ('–', '-', '—') else [],
        'relation': ('unassigned' if not labels else 'single_site' if len(labels) == 1 else
                     'overlapping' if len(values) > 1 else 'shared_shift'),
    }


def read_table(page, words, columns, expected):
    starts = [w['x0'] for _, w in columns]
    left = starts[0] - 2
    head = max(w['top'] for _, w in columns)
    last = columns[-1][1]
    table_rules = [line for line in page.lines if abs(line['top'] - line['bottom']) < 1 and
                   abs(line['x0'] - starts[0]) < 5 and line['x1'] >= last['x1'] - 2 and line['top'] > head]
    right = min([line['x1'] + .5 for line in table_rules if line['top'] - head < 15]
                or [min(page.width, last['x1'] + 12)])
    if expected:
        # The reviewed Wiley layout ends just before the adjacent prose column.
        right = min(right, 408)
    bottom = min([line['top'] for line in table_rules if line['top'] - head > 20] or [page.height])
    edges = [start - 2 for start in starts] + [right]
    coupling = next((word for key, word in columns if key == 'j'), None)
    j_hz = not coupling or 'hz' in normalize(coupling['text']) or any(
        abs(w['top'] - coupling['top']) < 3 and 0 <= w['x0'] - coupling['x1'] < 10
        and normalize(w['text']).strip('()[]') == 'hz' for w in words
    )
    shift_starts = [w['x0'] for key, w in columns if key in ('h', 'c')]
    anchors = sorted({round(w['top'], 2) for w in words if head + 5 < w['top'] < bottom
                      and any(abs(w['x0'] - start) < 2.6 for start in shift_starts)
                      and SHIFT.fullmatch(w['text'])})
    rows, issues = {'1H': [], '13C': []}, []
    heights = []
    for anchor in anchors:
        if heights and anchor - heights[-1] <= 2.2:
            continue
        if heights and anchor - heights[-1] > 24:
            break
        heights.append(anchor)
        if expected and len(heights) >= max(expected.values()):
            break
    for row, top in enumerate(heights, 1):
        chars = sorted((c for c in page.chars if top - 1.5 < c['top'] < top + 4.6
                        and left <= c['x0'] < right), key=lambda c: c['x0'])
        cells = {key: ''.join(c['text'] for c in chars if edges[i] <= c['x0'] < edges[i+1])
                 for i, (key, _) in enumerate(columns)}
        for shift, assign, nucleus in [('h','ha','1H'), ('c','ca','13C')]:
            if cells[shift]:
                try:
                    peak = parse_peak(cells[shift], cells[assign], nucleus, row, cells)
                    if peak['j_values_hz'] and not j_hz:
                        peak.update(j_values_reported=cells['j'], j_values_hz=[])
                    rows[nucleus].append(peak)
                except ValueError as exc:
                    issues.append(f'第 {row} 行 {nucleus} · {exc}')
            elif cells[assign]:
                issues.append(f'第 {row} 行 {nucleus} · 原子编号缺少位移')
    footnote_bottom = min(page.height, (heights[-1] if heights else head) + 92)
    notes = page.crop((left, (heights[-1] if heights else head) + 6, right, footnote_bottom)).extract_text() or ''
    for nucleus, peaks in rows.items():
        for mark in sorted({m for p in peaks for m in p['source_footnotes']}):
            members = [p for p in peaks if mark in p['source_footnotes']]
            pattern = rf'\[{mark}\]([^\[]*)'
            text = re.search(pattern, notes, re.I)
            if text and 'interchangeable' in text[1].casefold() and len(members) >= 2:
                for peak in members:
                    peak.update(relation='interchangeable', group=f'{nucleus}-{mark}')
            else:
                issues.append(f'{nucleus} 脚注 [{mark}] 需要复核')
        repeats = {}
        for peak in peaks:
            if len(peak['source_labels']) > 1 and not peak.get('group'):
                repeats.setdefault(tuple(peak['source_labels']), []).append(peak)
        for index, (labels, members) in enumerate(repeats.items(), 1):
            # Repeated lists may also be genuine overlap. Treat as a candidate set only.
            if len(members) == len(labels):
                for peak in members:
                    peak.update(relation='interchangeable', group=f'{nucleus}-repeat-{index}')
    return rows, issues, notes


def extract_document(path, profiles=None):
    import pdfplumber
    try:
        from .comparison_tables import ComparisonTables
    except ImportError:
        from comparison_tables import ComparisonTables

    profiles = profiles or {}
    tables, text_pages, doi = [], 0, None
    comparisons = ComparisonTables()
    with pdfplumber.open(path) as pdf:
        for number, page in enumerate(pdf.pages, 1):
            if not page.chars:
                continue
            text_pages += 1
            text = page.extract_text() or ''
            if re.search(r'natural\s*\(\s*δ', text, re.I) or comparisons.active:
                comparisons.read_page(page, number)
            if not doi and number <= 3:
                match = DOI.search(text)
                if match:
                    doi = match[0].rstrip('.,;)').lower()
            words = page.extract_words(x_tolerance=2.1, y_tolerance=2.2)
            for cols in headings(words):
                first = cols[0][1]
                caption = page.crop((max(0, first['x0'] - 7), max(0, first['top'] - 170),
                                     min(page.width, cols[-1][1]['x1'] + 12), first['top'] - 4)).extract_text() or ''
                matches = list(re.finditer(r'\bTable\s*(S?\d+)\s*[.:]?[^\n]*', caption, re.I))
                caption_text = matches[-1][0].strip() if matches else ''
                table = matches[-1][1].upper() if matches else None
                expected = profiles.get(str(number), {}).get('table_rows')
                rows, issues, notes = read_table(page, words, cols, expected)
                if expected and {k: len(v) for k, v in rows.items()} != expected:
                    issues.append('提取行数与已核对的原始表格不一致')
                if any(rows.values()) or issues:
                    tables.append({'page': number, 'table': table, 'caption': caption_text,
                                   'peaks': rows, 'issues': issues, 'footnotes_text': notes})
        return {'doi': doi, 'pages': len(pdf.pages), 'text_pages': text_pages,
                'tables': tables + list(comparisons.records.values())}


def probe_document(path):
    import pdfplumber

    with pdfplumber.open(path) as pdf:
        first = pdf.pages[0].extract_text() or ''
        second = (pdf.pages[1].extract_text() or '') if len(pdf.pages) > 1 else ''
        raw_title = str(pdf.metadata.get('Title') or '')
        opening = '\n'.join(first.splitlines()[:22])
        supplement = bool(re.search(r'\bsupporting information\b|\bsupplementary (?:information|material|data)\b',
                                    raw_title+'\n'+'\n'.join(first.splitlines()[:9]), re.I))
        doi_text = '\n'.join(str(pdf.metadata.get(k) or '') for k in ('DOI', 'Subject', 'Title')) + '\n' + opening
        if supplement:
            doi_text += '\n' + '\n'.join(second.splitlines()[:12])
        ids = {m[0].rstrip('.,;)').lower().split('#')[0] for m in DOI.finditer(doi_text)}
        doi = next(iter(ids)) if len(ids) == 1 else None
        title = re.sub(r'[^a-z0-9]+', ' ', raw_title.casefold()).strip()
        if re.search(r'^microsoft word|\.docx?$|supporting information|supplementary|^untitled', raw_title.casefold()):
            title = ''
        if supplement and not title:
            content = re.sub(r'(?i)\b(?:supporting information|supplementary (?:information|material|data))\b', '', opening)
            title = re.sub(r'[^a-z0-9]+', ' ', content.casefold()).strip()
        article = bool(re.search(r'\bAbstract\b|\bFull Paper\b|\bReceived\s*:', opening, re.I) or (doi and title))
        return {'doi': doi, 'title': title, 'supplement': supplement, 'article': article}


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf-8')
    path = Path(sys.argv[1])
    if sys.argv[2] == 'probe':
        print(json.dumps(probe_document(path), ensure_ascii=False))
        sys.exit()
    sample = json.loads(Path(__file__).with_name('mauritines.json').read_text(encoding='utf-8'))
    profiles = {str(c['pdf_page']): c for c in sample['compounds']} if sys.argv[2] == sample['pdf_sha256'] else {}
    print(json.dumps(extract_document(path, profiles), ensure_ascii=False))
