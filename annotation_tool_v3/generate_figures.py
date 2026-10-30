"""Generate editable SVG explanations from the implemented database and a real row."""
from html import escape
import json
import re
from pathlib import Path
from . import core as db

OUT = db.ROOT / "docs/database_v1_figures"
INK, MUTED, GREEN, BLUE, ORANGE = "#19382f", "#60716a", "#157d60", "#337a9c", "#af732e"


def text(x, y, value, size=22, color=INK, weight=400):
    return f'<text x="{x}" y="{y}" font-size="{size}" fill="{color}" font-weight="{weight}">{escape(str(value))}</text>'


def box(x, y, w, h, title, lines=(), accent=GREEN):
    s = f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="16" fill="white" stroke="#d4e2da" stroke-width="2"/><rect x="{x}" y="{y+18}" width="5" height="32" fill="{accent}"/>'
    s += text(x+24, y+43, title, 25, accent, 700)
    for i, line in enumerate(lines):
        s += text(x+24, y+83+i*32, line, 20)
    return s


def arrow(points, color=GREEN):
    return f'<polyline points="{points}" fill="none" stroke="{color}" stroke-width="3" marker-end="url(#arrow)"/>'


def svg(name, title, subtitle, body, height=860):
    content = f'''<svg xmlns="http://www.w3.org/2000/svg" width="1440" height="{height}" viewBox="0 0 1440 {height}" role="img" aria-label="{escape(title)}">
<defs><marker id="arrow" markerWidth="9" markerHeight="9" refX="8" refY="4" orient="auto"><path d="M0,0 L8,4 L0,8" fill="{GREEN}"/></marker></defs>
<rect width="1440" height="{height}" fill="#f4f8f4"/>
<g font-family="Microsoft YaHei, Noto Sans CJK SC, Arial, sans-serif">{text(52,68,title,34,INK,700)}{text(52,108,subtitle,21,MUTED)}{body}
{text(52,height-26,'实验 NMR 归属数据库 · v1 / 2026-09-29 · 矢量图，可缩放与编辑',17,MUTED)}</g></svg>'''
    (OUT/name).write_text(content, encoding="utf-8")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    body = box(52,155,410,165,"NMRShiftDB2 · 实验谱",["按每条谱筛除 HOSE / ACD 等预测", "保留原始归属、原子映射和全部条件", "同一分子的多次实验分别保存"])
    body += box(52,355,410,165,"NMRPeak · 去重目录",["只作为 NMRexp 的结构候选入口", "与实验 ShiftDB2 按完整结构去重", "不把 LMDB 序号当原始数据行号"],BLUE)
    body += box(530,355,400,165,"回到 NMRexp 原始表",["结构、峰指纹、溶剂和文献核对", "恢复 DOI、原文、原表行与页码", "定位不唯一：留在迁移目录，不猜"],BLUE)
    body += box(1000,355,388,165,"NMRexp · 实验谱",["重复但缺谱：条件相容才补谱", "新结构：H / C 双谱待标注", "不把既有标签移植到另一条实验"],BLUE)
    body += box(530,155,858,165,"SQLite · 日常工作的主数据库",["一个分子 → 多条实验谱 → 每条谱的信号 → 每个信号的归属", "保留原库标签；新增人工标注。没有复核、DFT、预测或 FID 模块。", "每次保存均有原子编号校验和并发版本检查"])
    body += arrow("462,230 520,230") + arrow("462,440 520,440") + arrow("930,440 990,440") + arrow("1195,355 1195,330")
    body += box(52,580,420,165,"工作顺序",["① 已有部分归属 → 补剩余谱 / 峰", "② 两谱均未归属 → 逐步新标", "缺实验谱单列；完整记录保留"])
    body += box(525,580,405,165,"随时导出",["JSONL：完整、可读、可交换", "Parquet：列式分析与批量处理", "LMDB：训练消费，带归属掩码"])
    body += box(983,580,405,165,"不做隐含假设",["溶剂相同 ≠ 温度或样品相同", "H/C 的 MHz 数字不应直接比相等", "缺失条件是 null，不是默认值"],ORANGE)
    body += arrow("960,320 960,555 745,555 745,568")
    svg("01_build_flow.svg","从两个数据来源，到一套日常标注数据库","原始文件只读；新库独立保存。完整结构去重不合并立体异构体、盐形式或不同电荷。",body)

    body = box(52,160,410,215,"molecules · 分子",["主键 id；canonical SMILES", "InChIKey / 分子式 / 分子量", "固定的原子编号与结构对称组", "任务队列 / 备注 / 保存版本"])
    body += box(530,160,410,215,"spectra · 实验谱",["主键 id；外键 molecule_id", "核种 / 来源 / 各自的实验条件", "peaks_json：信号与原始数值", "source_json：原表完整数据"])
    body += box(1008,160,380,215,"annotations · 归属",["联合键 spectrum_id + peak_id", "状态 / 原子编号集合 / 备注", "原库或人工 / 保存时间", "标签属于具体实验，不属于分子"])
    body += arrow("462,256 520,256")+text(471,235,"1:N",18,MUTED)
    body += arrow("940,256 998,256")+text(949,235,"1:N",18,MUTED)
    body += box(530,450,410,190,"pairs · H/C 配对",["同一个 molecule_id", "c_id + h_id 指向两条具体谱", "配对依据 / 溶剂一致性", "不知道样品身份时明确写未知"])
    body += arrow("735,375 735,438")
    body += box(52,450,410,255,"迁移辅助表 · 非标注流程",["candidates：处理后目录与回源状态", "source_rows：原始 NMRexp 行", "progress：中断后继续的位置", "meta：输入文件与校验信息", "legacy_annotations：旧文件原样留存", "import_issues：程序导入异常"],BLUE)
    body += box(1008,450,380,255,"为什么不是一张大表？",["一个分子可能有多种溶剂", "同一溶剂也可能有多次实验", "C/H 归属进度可以不同", "一个信号可以有多个归属原子", "歧义不被强制压成一个答案", "出口可再聚合为每分子一条"])
    body += text(530,703,"日常主表只有：分子、实验谱、归属、配对。",23,GREEN,700)
    svg("02_database_structure.svg","数据库的基本组成单元","SQLite 管理关系；JSON 字段保存可扩展的峰信息与原始列。新增信息不需要改写历史原文。",body)

    c = db.connect()
    row = c.execute("SELECT id FROM molecules WHERE smiles='CCO'").fetchone()
    r = db.record(c,row[0])
    c.close()
    (OUT/"example_ethanol.json").write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding="utf-8")
    chosen = {s["nucleus"]:s for s in r["spectra"] if s["solvent"]=="D2O"}
    body = box(52,160,480,500,"分子层 · 乙醇 C₂H₆O",["SMILES：CCO", "原子编号：C0 — C1 — O2", "H 标签记录在所连接的重原子上"])
    drawing = db.molecule_svg(r["smiles"])
    drawing = drawing[drawing.index("<svg"):]
    drawing = re.sub(r"<svg[^>]*>",'<svg x="64" y="315" width="456" height="278" viewBox="0 0 720 440">',drawing,count=1)
    body += drawing
    body += text(76,618,"这一分子在库中实际保存 5 条实验谱。",20,MUTED)
    for nuc,x,title in (("13C",585,"¹³C · D₂O"),("1H",1000,"¹H · D₂O")):
        s=chosen[nuc]
        body += box(x,160,388,365,title,["来源：NMRShiftDB2 10009222", "核种条件独立保存；下列为原库标签"])
        body += text(x+24,310,"δ / ppm",21,MUTED,700)+text(x+205,310,"归属目标",21,MUTED,700)
        for i,p in enumerate(s["peaks"]):
            targets=" / ".join(r["structure"]["atoms"][a]["symbol"]+str(a) for a in p["assignment"]["atom_indices"])
            body += text(x+24,361+57*i,p["shift_ppm"],30,INK,700)+text(x+205,361+57*i,targets,28,GREEN,700)
        body += text(x+24,496,"归属来源：原库导入，非本次人工填写",18,MUTED)
    body += box(585,565,803,168,"这条记录怎样理解？",["D₂O 中的 C0：δC 19.49 ppm；连接 C0 的 H：δH 1.17 ppm。", "C1：δC 60.13 ppm；连接 C1 的 H：δH 3.641 ppm。", "没有列出的 OH 信号不填 0、不编造归属；同一分子另有其他条件的谱。"])
    body += text(52,785,"图中只展示 D₂O 配对；完整示例 JSON 保留全部 5 条谱、来源字段和固定原子编号。",22,MUTED)
    svg("03_real_molecule_example.svg","一个真实分子，怎样落在数据库里","以下位移直接读取本地已迁移数据，不是示意数值。结构图及表格均为矢量。",body,880)
    print(str(OUT))


if __name__ == "__main__":
    main()
