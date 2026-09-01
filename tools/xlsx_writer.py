"""极简 xlsx 写出器（纯标准库实现，不依赖 openpyxl / xlsxwriter）。

项目运行环境未安装第三方 Excel 库，且构建机可能离线，因此这里直接按 OOXML
最小可用集合拼装 .xlsx：内容类型清单、工作簿、若干工作表、样式表。字符串统一
使用 inline string（省去共享字符串表），数值直接写 <v>。

对外只暴露 Workbook.add_sheet / Workbook.save 两个入口，样式用常量 STYLE_* 指定。
被 tools/test_report.py 用于输出“按检查类别分 sheet”的测试报表。
"""

import re
import zipfile

# 工作表单元格样式（索引与 styles.xml 的 cellXfs 顺序一一对应）
STYLE_DEFAULT = 0
STYLE_HEADER = 1
STYLE_ERROR = 2
STYLE_WARN = 3
STYLE_INFO = 4
STYLE_TITLE = 5
STYLE_WRAP = 6
STYLE_OK = 7

_ESCAPES = {"&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&apos;"}


def _esc(text):
    """XML 转义，并剔除 XML 1.0 不允许的控制字符。"""
    out = []
    for ch in str(text):
        if ch in _ESCAPES:
            out.append(_ESCAPES[ch])
        elif ch in ("\t", "\n", "\r"):
            out.append(" ")
        elif ord(ch) < 0x20:
            continue
        else:
            out.append(ch)
    return "".join(out)


def _col_name(index):
    """0 -> A, 25 -> Z, 26 -> AA。"""
    name = ""
    index += 1
    while index:
        index, rem = divmod(index - 1, 26)
        name = chr(ord("A") + rem) + name
    return name


def sanitize_sheet_name(name):
    """Excel 工作表名限制：禁用 []:*?/\\，长度不超过 31 字符。"""
    cleaned = re.sub(r"[\[\]:*?/\\]", "_", str(name)).strip() or "sheet"
    return cleaned[:31]


class _Sheet(object):
    def __init__(self, name, headers, rows, widths, row_styles, freeze, autofilter):
        self.name = name
        self.headers = list(headers)
        self.rows = [list(r) for r in rows]
        self.widths = list(widths) if widths else None
        self.row_styles = list(row_styles) if row_styles else None
        self.freeze = freeze
        self.autofilter = autofilter


class Workbook(object):
    """按 sheet 顺序累积数据，最后一次性写出 .xlsx。"""

    def __init__(self):
        self._sheets = []

    def add_sheet(self, name, headers, rows, widths=None, row_styles=None,
                  freeze=True, autofilter=True):
        """追加一张工作表。

        name        工作表名（自动裁剪到 31 字符并去重）
        headers     表头文本列表；空列表表示无表头
        rows        每行一个列表，元素为 str / int / float / None
        widths      列宽（字符数）列表，可为 None
        row_styles  与 rows 等长的样式常量列表，用于按严重度着色
        freeze      是否冻结表头行
        autofilter  是否为表头加筛选器
        """
        final = sanitize_sheet_name(name)
        existing = set(s.name for s in self._sheets)
        if final in existing:
            base = final[:28]
            for i in range(2, 100):
                cand = "%s_%d" % (base, i)
                if cand not in existing:
                    final = cand
                    break
        self._sheets.append(
            _Sheet(final, headers, rows, widths, row_styles, freeze, autofilter))
        return final

    def save(self, path):
        """写出 .xlsx；无 sheet 时补一张空表，避免生成 Excel 打不开的空工作簿。"""
        if not self._sheets:
            self.add_sheet("empty", ["(no data)"], [])
        with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as zf:
            zf.writestr("[Content_Types].xml", self._content_types())
            zf.writestr("_rels/.rels", _ROOT_RELS)
            zf.writestr("xl/workbook.xml", self._workbook_xml())
            zf.writestr("xl/_rels/workbook.xml.rels", self._workbook_rels())
            zf.writestr("xl/styles.xml", _STYLES_XML)
            for i, sheet in enumerate(self._sheets):
                zf.writestr("xl/worksheets/sheet%d.xml" % (i + 1),
                            self._sheet_xml(sheet))

    def _content_types(self):
        parts = [
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
            '<Types xmlns="http://schemas.openxmlformats.org/package/2006/'
            'content-types">',
            '<Default Extension="rels" ContentType="application/'
            'vnd.openxmlformats-package.relationships+xml"/>',
            '<Default Extension="xml" ContentType="application/xml"/>',
            '<Override PartName="/xl/workbook.xml" ContentType="application/'
            'vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>',
            '<Override PartName="/xl/styles.xml" ContentType="application/'
            'vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>',
        ]
        for i in range(len(self._sheets)):
            parts.append(
                '<Override PartName="/xl/worksheets/sheet%d.xml" '
                'ContentType="application/vnd.openxmlformats-officedocument.'
                'spreadsheetml.worksheet+xml"/>' % (i + 1))
        parts.append("</Types>")
        return "".join(parts)

    def _workbook_xml(self):
        sheets = "".join(
            '<sheet name="%s" sheetId="%d" r:id="rId%d"/>' % (_esc(s.name), i + 1, i + 1)
            for i, s in enumerate(self._sheets))
        return ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                '<workbook xmlns="http://schemas.openxmlformats.org/'
                'spreadsheetml/2006/main" xmlns:r="http://schemas.'
                'openxmlformats.org/officeDocument/2006/relationships">'
                '<sheets>%s</sheets></workbook>' % sheets)

    def _workbook_rels(self):
        rels = "".join(
            '<Relationship Id="rId%d" Type="http://schemas.openxmlformats.org/'
            'officeDocument/2006/relationships/worksheet" '
            'Target="worksheets/sheet%d.xml"/>' % (i + 1, i + 1)
            for i in range(len(self._sheets)))
        rels += ('<Relationship Id="rId%d" Type="http://schemas.openxmlformats.'
                 'org/officeDocument/2006/relationships/styles" '
                 'Target="styles.xml"/>' % (len(self._sheets) + 1))
        return ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                '<Relationships xmlns="http://schemas.openxmlformats.org/'
                'package/2006/relationships">%s</Relationships>' % rels)

    def _sheet_xml(self, sheet):
        ncols = max([len(sheet.headers)] + [len(r) for r in sheet.rows] + [1])
        parts = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                 '<worksheet xmlns="http://schemas.openxmlformats.org/'
                 'spreadsheetml/2006/main">']
        if sheet.freeze and sheet.headers:
            parts.append('<sheetViews><sheetView workbookViewId="0">'
                         '<pane ySplit="1" topLeftCell="A2" activePane="bottomLeft"'
                         ' state="frozen"/></sheetView></sheetViews>')
        if sheet.widths:
            cols = "".join(
                '<col min="%d" max="%d" width="%.1f" customWidth="1"/>'
                % (i + 1, i + 1, float(w))
                for i, w in enumerate(sheet.widths) if w)
            if cols:
                parts.append("<cols>%s</cols>" % cols)
        parts.append("<sheetData>")

        row_no = 1
        if sheet.headers:
            parts.append(self._row_xml(row_no, sheet.headers, STYLE_HEADER))
            row_no += 1
        for i, row in enumerate(sheet.rows):
            style = sheet.row_styles[i] if sheet.row_styles else STYLE_DEFAULT
            parts.append(self._row_xml(row_no, row, style))
            row_no += 1
        parts.append("</sheetData>")

        if sheet.autofilter and sheet.headers and sheet.rows:
            parts.append('<autoFilter ref="A1:%s%d"/>'
                         % (_col_name(ncols - 1), len(sheet.rows) + 1))
        parts.append("</worksheet>")
        return "".join(parts)

    @staticmethod
    def _row_xml(row_no, values, style):
        cells = []
        for col, value in enumerate(values):
            if value is None or value == "":
                continue
            ref = "%s%d" % (_col_name(col), row_no)
            if isinstance(value, bool):
                value = "TRUE" if value else "FALSE"
            if isinstance(value, (int, float)):
                num = round(float(value), 4)
                if isinstance(value, int):
                    num = value
                cells.append('<c r="%s" s="%d"><v>%s</v></c>'
                             % (ref, style, repr(num) if isinstance(num, float)
                                else str(num)))
            else:
                cells.append('<c r="%s" s="%d" t="inlineStr"><is><t xml:space='
                             '"preserve">%s</t></is></c>' % (ref, style, _esc(value)))
        return '<row r="%d">%s</row>' % (row_no, "".join(cells))


_ROOT_RELS = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
    '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/'
    'relationships"><Relationship Id="rId1" Type="http://schemas.'
    'openxmlformats.org/officeDocument/2006/relationships/officeDocument" '
    'Target="xl/workbook.xml"/></Relationships>')

# 字体/填充/边框/单元格格式的顺序即 STYLE_* 常量的取值来源。
_STYLES_XML = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
    '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
    '<fonts count="7">'
    '<font><sz val="11"/><name val="Calibri"/></font>'
    '<font><b/><sz val="11"/><color rgb="FFFFFFFF"/><name val="Calibri"/></font>'
    '<font><b/><sz val="11"/><color rgb="FFC00000"/><name val="Calibri"/></font>'
    '<font><sz val="11"/><color rgb="FFBF6000"/><name val="Calibri"/></font>'
    '<font><sz val="11"/><color rgb="FF7F7F7F"/><name val="Calibri"/></font>'
    '<font><b/><sz val="13"/><name val="Calibri"/></font>'
    '<font><sz val="11"/><color rgb="FF00703C"/><name val="Calibri"/></font>'
    '</fonts>'
    '<fills count="3">'
    '<fill><patternFill patternType="none"/></fill>'
    '<fill><patternFill patternType="gray125"/></fill>'
    '<fill><patternFill patternType="solid"><fgColor rgb="FF4472C4"/>'
    '<bgColor indexed="64"/></patternFill></fill>'
    '</fills>'
    '<borders count="2"><border/>'
    '<border><bottom style="thin"><color rgb="FF9E9E9E"/></bottom></border>'
    '</borders>'
    '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/>'
    '</cellStyleXfs>'
    '<cellXfs count="8">'
    '<xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>'
    '<xf numFmtId="0" fontId="1" fillId="2" borderId="1" xfId="0" applyFont="1"'
    ' applyFill="1" applyBorder="1"><alignment vertical="center"/></xf>'
    '<xf numFmtId="0" fontId="2" fillId="0" borderId="0" xfId="0" applyFont="1"/>'
    '<xf numFmtId="0" fontId="3" fillId="0" borderId="0" xfId="0" applyFont="1"/>'
    '<xf numFmtId="0" fontId="4" fillId="0" borderId="0" xfId="0" applyFont="1"/>'
    '<xf numFmtId="0" fontId="5" fillId="0" borderId="0" xfId="0" applyFont="1"/>'
    '<xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0" applyAlignment="1">'
    '<alignment wrapText="1" vertical="top"/></xf>'
    '<xf numFmtId="0" fontId="6" fillId="0" borderId="0" xfId="0" applyFont="1"/>'
    '</cellXfs>'
    '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/>'
    '</cellStyles></styleSheet>')


SEVERITY_STYLE = {
    "error": STYLE_ERROR,
    "warn": STYLE_WARN,
    "info": STYLE_INFO,
    "pass": STYLE_OK,
}


def style_for_severity(severity):
    """严重度 -> 单元格样式，未知severity按默认样式处理。"""
    return SEVERITY_STYLE.get(str(severity).lower(), STYLE_DEFAULT)
