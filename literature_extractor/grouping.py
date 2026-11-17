"""Identify and move loose PDFs into one folder per confidently matched article."""

from difflib import SequenceMatcher
from pathlib import Path
import os
import re

from .extract import COMPARISON_SAMPLE, SAMPLE, digest, inspect_pdf, source_paths

SUPPLEMENT = re.compile(r'support(?:ing|ingmaterials)|supplement(?:ary|al)?|(?:^|[^a-z0-9])(?:si|esm|s\d{1,2})(?:$|[^a-z0-9])', re.I)
KNOWN = {p[role + '_sha256']: (p['doi'], role) for p in (SAMPLE, COMPARISON_SAMPLE)
         for role in ('pdf', 'supplement')}


def stem(path):
    cleaned = SUPPLEMENT.sub(' ', Path(path).stem.casefold())
    return re.sub(r'[^a-z0-9]+', ' ', cleaned).strip()


def score(supp, paper):
    if supp['doi'] and paper['doi']:
        return (4, 'DOI 一致') if supp['doi'] == paper['doi'] else (0, '')
    a, b = supp['title'], paper['title']
    if min(len(a), len(b)) >= 35 and SequenceMatcher(None, a, b, autojunk=False).ratio() >= .92:
        return 3, '首页题名一致'
    a, b = stem(supp['file']), stem(paper['file'])
    if len(a) >= 10 and a == b:
        return 2, '文件主名称一致'
    return 0, ''


def identify(path, code):
    if code in KNOWN:
        doi, role = KNOWN[code]
        return {'file': path.name, 'sha256': code, 'doi': doi, 'title': '', 'role': role}
    try:
        metadata = inspect_pdf(path)
        supplement = metadata['supplement'] or SUPPLEMENT.search(path.stem) or re.fullmatch(r'[a-z0-9]+_s', path.stem, re.I)
        role = 'supplement' if supplement else 'pdf' if metadata['article'] else 'unknown'
        return {'file': path.name, 'sha256': code, 'doi': metadata['doi'], 'title': metadata['title'], 'role': role}
    except Exception as exc:
        return {'file': path.name, 'sha256': code, 'role': 'unknown', 'reason': str(exc)}


def group_target(folder, paper):
    name = re.sub(r'[^\w.-]+', '_', paper['doi'] or Path(paper['file']).stem, flags=re.UNICODE)[:72].strip('._ ')
    return folder / f'{name or "article"}_{paper["sha256"][:8]}'


def relocate(source, destination):
    if os.name == 'nt':
        source.rename(destination)  # Windows rejects an existing destination.
    else:
        os.link(source, destination)  # Publish without overwriting an existing file.
        source.unlink()


def grouping_plan(folder, *, apply=False):
    folder = folder.resolve()
    loose = [path for path in source_paths(folder) if path.suffix.lower() == '.pdf']
    if not loose:
        return {'groups': [], 'pending': [], 'moved': 0, 'duplicates': 0}
    hashes, metadata = {}, {}

    def checksum(path):
        if path not in hashes:
            hashes[path] = digest(path)
        return hashes[path]

    def identify_once(path):
        code = checksum(path)
        if code not in metadata or metadata[code]['role'] == 'unknown':
            metadata[code] = identify(path, code)
        return {**metadata[code], 'file': path.name, 'path': path}

    papers, groups, pending = {}, {}, []
    for child in sorted(folder.iterdir()):
        if (not re.search(r'_[a-f0-9]{8}$', child.name) or child.is_symlink()
                or not child.is_dir() or child.resolve().parent != folder):
            continue
        for path in source_paths(child):
            if path.suffix.lower() == '.pdf' and checksum(path)[:8] == child.name[-8:]:
                paper = identify_once(path)
                if paper['role'] == 'pdf':
                    papers.setdefault(paper['sha256'], (paper, child))

    attachments = []
    for path in loose:
        doc = identify_once(path)
        key = doc['sha256']
        if doc['role'] == 'pdf' or key in papers:
            papers.setdefault(key, (doc, group_target(folder, doc)))
            groups.setdefault(key, []).append((doc, '正文'))
        elif doc['role'] == 'supplement':
            attachments.append(doc)
        else:
            pending.append({'file': doc['file'], 'reason': doc.get('reason') or '无法可靠区分正文和附件'})

    for doc in attachments:
        matches = [(score(doc, paper), key) for key, (paper, _) in papers.items()]
        best = max((points for (points, _), _ in matches), default=0)
        winners = [(basis, key) for (points, basis), key in matches if points == best and points]
        if len(winners) == 1:
            basis, key = winners[0]
            groups.setdefault(key, []).append((doc, basis))
        else:
            pending.append({'file': doc['file'], 'reason': '匹配到多个候选正文' if winners else '没有可靠的正文匹配'})

    results, moved, duplicates = [], 0, 0
    for key, entries in groups.items():
        paper, target = papers[key]
        result = {'folder': target.name, 'paper': paper['file'], 'files': [doc['file'] for doc, _ in entries],
                  'evidence': list(dict.fromkeys(basis for _, basis in entries))}
        results.append(result)
        result['status'] = '可分组'
        if target.is_symlink() or (target.exists() and (not target.is_dir() or target.resolve().parent != folder)):
            result['status'] = '目标目录存在冲突，请检查'
            continue
        existing = {checksum(path): path for path in source_paths(target) if path.suffix.lower() == '.pdf'} if target.exists() else {}
        if not apply:
            continue
        target.mkdir(exist_ok=True)
        for doc, _ in entries:
            source, code = doc['path'], doc['sha256']
            destination = target / source.name
            try:
                if source.is_symlink() or digest(source) != code:
                    raise ValueError('分组过程中源文件发生变更')
                if code in existing and not existing[code].is_symlink() and digest(existing[code]) == code:
                    source.unlink()
                    duplicates += 1
                else:
                    relocate(source, destination)
                    existing[code] = destination
                    moved += 1
            except (OSError, ValueError) as exc:
                pending.append({'file': doc['file'], 'reason': f'移动失败 · {exc}'})
                result['status'] = '部分文件待确认'
        if result['status'] == '可分组':
            result['status'] = '已分组'
    return {'groups': results, 'pending': pending, 'moved': moved, 'duplicates': duplicates}
