"""Draw actual stored prediction examples as vector structures and PNG previews."""
from pathlib import Path
import json
from html import escape

from PIL import Image,ImageDraw,ImageFont
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/atompeak_partial_scratch_v1/explanation_20261031'


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    rows={r['id']:r for line in (ROOT/'results/atompeak_partial_scratch_v1/test_best_generation_v1.jsonl').read_text(encoding='utf-8').splitlines() if (r:=json.loads(line))}
    ranking=rows['M21f87339c9c7e755eefbd26e'];isomer=rows['M2676ddbf02b0003a31baff71']
    W,H=1500,1110
    canvas=Image.new('RGB',(W,H),'#f4f7f8');pen=ImageDraw.Draw(canvas)
    svg=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
         f'<rect width="{W}" height="{H}" fill="#f4f7f8"/>']
    font_path='C:/Windows/Fonts/msyh.ttc'
    def text(x,y,s,size=24,color='#243840'):
        pen.text((x,y),s,font=ImageFont.truetype(font_path,size),fill=color)
        svg.append(f'<text x="{x}" y="{y+size}" fill="{color}" font-family="Microsoft YaHei,sans-serif" font-size="{size}">{escape(s)}</text>')
    def card(x,y,title,smiles,footer,accent='#176d69'):
        width,height=460,330
        pen.rounded_rectangle((x,y,x+width,y+height),radius=14,fill='white',outline='#d5dfe2',width=2)
        svg.append(f'<rect x="{x}" y="{y}" width="{width}" height="{height}" rx="14" fill="white" stroke="#d5dfe2"/>')
        text(x+20,y+14,title,25,accent)
        mol=Chem.MolFromSmiles(smiles)
        drawer=rdMolDraw2D.MolDraw2DSVG(420,180);drawer.DrawMolecule(mol);drawer.FinishDrawing()
        tile=drawer.GetDrawingText();tile=tile[tile.index('<svg'):]
        svg.append(f'<g transform="translate({x+20},{y+57})">{tile}</g>')
        raster=rdMolDraw2D.MolDraw2DCairo(420,180);raster.DrawMolecule(mol);raster.FinishDrawing()
        import io
        canvas.paste(Image.open(io.BytesIO(raster.GetDrawingText())).convert('RGB'),(x+20,y+57))
        text(x+20,y+245,rdMolDescriptors.CalcMolFormula(mol),29,accent)
        text(x+20,y+290,footer,20)
    text(30,18,'实际预测示例：能画出结构 ≠ 分子式符合 ≠ 找对连接关系',32)
    text(30,65,'A. 已知输入 C3H5Cl；正确答案在第二位，首选多了一个 Cl、少了一个 H。',24)
    card(30,110,'真实目标',ranking['target_for_evaluation_only'],'M21f87339…')
    card(520,110,'原候选第 1 位：分子式错',ranking['beams'][0]['smiles'],'可解析，但不符合已知组成','#b24f3c')
    card(1010,110,'原候选第 2 位：正确',ranking['beams'][1]['smiles'],'分子式检查可使其升到第一')
    text(30,466,'B. 已知输入 C3H9N；下面两个都符合分子式，但原子的连接方式不同。',24)
    card(30,510,'真实目标：乙基甲基胺',isomer['target_for_evaluation_only'],'正确结构原本排在第 4 位')
    card(520,510,'原候选第 1 位：丙胺',isomer['beams'][0]['smiles'],'分子式检查无法区分这两者','#b24f3c')
    x,y=1010,510
    pen.rounded_rectangle((x,y,x+460,y+330),radius=14,fill='white',outline='#d5dfe2',width=2)
    svg.append(f'<rect x="{x}" y="{y}" width="460" height="330" rx="14" fill="white" stroke="#d5dfe2"/>')
    text(x+20,y+14,'无法解析的实际候选',25,'#b24f3c')
    text(x+20,y+90,'Cl=C/C=C/C',34,'#b24f3c')
    text(x+20,y+158,'Cl 被写成双键连接；解析器拒绝。',22)
    text(x+20,y+201,'不能当作有效分子画出结构式，',22)
    text(x+20,y+240,'也不能拿它参与化学结构比较。',22)
    text(30,880,'C. “不同字符串、相同结构”的教学例子：CCO 与 OCC 都是乙醇，只算一个候选。',24)
    text(30,927,'上方 A 的实际记录：请求 10 个 → 保存 9 个 → 7 个合法不同结构 + 2 个非法字符串。',24)
    text(30,974,'按现有生成逻辑推算，另 1 个重复结构已被删除；旧输出没有保存被删的原字符串。',23)
    text(30,1032,'绿色表示目标/正确候选；红色表示错误。图片来自这批预测，不是重新生成模型结果。',22)
    svg.append('</svg>')
    (OUT/'candidate_examples.svg').write_text('\n'.join(svg),encoding='utf-8')
    canvas.save(OUT/'candidate_examples.png')
    print(OUT/'candidate_examples.svg')


if __name__=='__main__':main()
