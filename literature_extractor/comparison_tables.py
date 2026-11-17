"""NMR tables with position, natural sample and synthetic sample columns."""

import re


SITE = re.compile(r'([1-9]\d*)(?:\s*-\s*(NH))?$', re.I)
DECIMAL = r'\d+(?:\.\d+)'
PROTON_CELL = re.compile(
    rf'(?P<shifts>{DECIMAL}(?:\s*,\s*{DECIMAL})*)\s*,?\s*(?:\((?P<detail>[^()]*)\)|,\s*(?P<bare>[a-z]+))?', re.I)


def site(raw):
    return SITE.fullmatch((raw or '').strip())


def continuation_rows(page, grid):
    """Recover text above the first complete rectangle of an open-top grid."""
    left, top, _, _ = grid.bbox
    starts = [edge['top'] for edge in page.edges if edge['orientation'] == 'v' and
              abs(edge['x0'] - left) < 2 and edge['top'] < top - 4 < edge['bottom']]
    if not starts or not 4 < top - max(starts) < 38:
        return []
    bounds = grid.rows[0].cells
    if any(cell is None for cell in bounds):
        return []
    cells = [page.crop((x0, max(starts) - .25, x1, top)).extract_text() or '' for x0, _, x1, _ in bounds]
    return [cells] if site(cells[0]) else []


def detect_header(rows):
    if len(rows) < 2 or len(rows[0]) != 4 or len(rows[1]) != 4:
        return None
    if (rows[0][0] or '').strip().lower() not in ('no.', 'position', 'atom'):
        return None
    headings = [(text or '').lower() for text in rows[1]]
    if 'natural' not in headings[1] or 'synthetic' not in headings[2] or 'δ' not in headings[2]:
        return None
    name = (rows[0][1] or '').strip()
    if not name:
        return None
    if re.search(r'δ\s*h|\b1h\b', headings[2]):
        return name, '1H'
    populated = [(row[2] or '').strip() for row in rows[2:] if site(row[0]) and row[2]]
    if re.search(r'δ\s*c|\b13c\b', headings[2]) or populated and all(re.fullmatch(DECIMAL, cell) for cell in populated):
        return name, '13C'
    return None


def read_cells(raw, nucleus):
    if nucleus == '13C':
        if not re.fullmatch(DECIMAL, raw):
            raise ValueError('碳位移格式无法识别')
        return [(raw, None, [])]
    peaks = []
    for line in raw.splitlines():
        cell = PROTON_CELL.fullmatch(line.strip())
        if not cell:
            raise ValueError('氢位移、峰型或 J 值格式无法识别')
        detail = (cell['detail'] or cell['bare'] or '').split(',')
        multiplicity = detail[0].strip() or None
        j_values = [float(part.strip()) for part in detail[1:]]
        for delta in cell['shifts'].split(','):
            peaks.append((delta.strip(), multiplicity, j_values))
    return peaks


class ComparisonTables:
    def __init__(self):
        self.active = None
        self.records = {}

    def read_page(self, page, page_number):
        for grid in page.find_tables():
            rows = grid.extract(x_tolerance=1)
            header = detect_header(rows)
            if header:
                name, nucleus = header
                record = self.records.setdefault(name, {
                    'page': page_number, 'table': None, 'caption': name, 'comparison': 'synthetic',
                    'component_pages': {}, 'peaks': {'1H': [], '13C': []}, 'issues': [], 'footnotes_text': '',
                })
                self.active = (record, nucleus, page_number, grid.bbox[0], 0)
                self.read_grid(record, nucleus, page_number, rows[2:])
            elif self.active and rows and len(rows[0]) == 4 and site(rows[0][0]):
                record, nucleus, previous, left, last = self.active
                if page_number == previous + 1 and abs(grid.bbox[0] - left) < 3 and grid.bbox[1] < 150:
                    data = continuation_rows(page, grid) + rows
                    if int(site(data[0][0])[1]) != last + 1:
                        record['issues'].append(f'第 {page_number} 页的续表位点编号不连续')
                    self.read_grid(record, nucleus, page_number, data)

    def read_grid(self, record, nucleus, page_number, rows):
        peaks = record['peaks'][nucleus]
        record['component_pages'].setdefault(nucleus, []).append(page_number)
        last = self.active[4]
        for columns in rows:
            if len(columns) != 4 or not site(columns[0]):
                continue
            atom, nh = site(columns[0]).groups()
            last = int(atom)
            raw = (columns[2] or '').strip()
            if not raw:
                continue
            label = f"{'H' if nucleus == '1H' else 'C'}{atom}{'-NH' if nh else ''}"
            try:
                for delta, mult, j_values in read_cells(raw, nucleus):
                    peaks.append({'row': len(peaks) + 1, 'shift_text': delta, 'shift_ppm': float(delta),
                                  'range_ppm': None, 'source_labels': [label], 'source_assignment_text': label,
                                  'source_footnotes': [], 'multiplicity': mult, 'reported_nH': None,
                                  'j_values_hz': j_values, 'relation': 'single_site', 'source_page': page_number,
                                  'source_row_label': columns[0], 'synthetic_cell': raw,
                                  'natural_reference_cell': columns[1], 'difference_cell': columns[3]})
            except ValueError as exc:
                record['issues'].append(f'第 {page_number} 页 {label} · {exc}（{raw}）')
        self.active = (record, nucleus, page_number, self.active[3], last)
