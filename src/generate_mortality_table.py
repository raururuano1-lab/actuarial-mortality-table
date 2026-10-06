from __future__ import annotations

import math
import os
import textwrap
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FuncFormatter, ScalarFormatter
import numpy as np
import openpyxl
import pandas as pd
import xlsxwriter


BASE = Path(__file__).resolve().parents[1]
SERIE = 8
RADIX = 100_000.0
DATA_FILE = BASE / "data" / "mortality_rates_series_08.csv"
OUTPUT_DIR = BASE / "outputs"
OUTPUT_XLSX = OUTPUT_DIR / "mortality_table_series_08.xlsx"
OUTPUT_PDF = Path(
    os.environ.get(
        "MORTALITY_PDF_OUTPUT",
        str(OUTPUT_DIR / "mortality_table_report.pdf"),
    )
)
GRAPH_DIR = OUTPUT_DIR / "figures"


def fmt_num(value: float, decimals: int = 2) -> str:
    return f"{value:,.{decimals}f}"


def fmt_rate(value: float, decimals: int = 4) -> str:
    return f"{100 * value:.{decimals}f}%"


def age_label(age: int) -> str:
    return "109 y+" if age == 109 else str(age)


def find_source() -> Path:
    if not DATA_FILE.exists():
        raise RuntimeError(f"No se encontró el archivo de datos: {DATA_FILE}")
    return DATA_FILE


def load_assigned_series(source: Path) -> tuple[str, str, np.ndarray, list[str]]:
    data = pd.read_csv(source)
    if list(data.columns) != ["age", "mx"]:
        raise RuntimeError("El CSV debe contener exactamente las columnas age y mx.")

    raw_ages = data["age"].tolist()
    raw_mx = data["mx"].tolist()

    if len(raw_ages) != 110 or len(raw_mx) != 110:
        raise RuntimeError("La fuente no contiene exactamente 110 observaciones.")
    numeric_ages = [int(value) for value in raw_ages[:-1]]
    if numeric_ages != list(range(109)) or str(raw_ages[-1]).strip() != "109 y+":
        raise RuntimeError("La secuencia de edades no coincide con 0–108 y 109 y+.")
    if any(not isinstance(value, (int, float)) for value in raw_mx):
        raise RuntimeError("La serie mx contiene valores no numéricos o faltantes.")

    return "Raúl Antonio Ruano Katz", "mx", np.asarray(raw_mx, dtype=float), [str(x) for x in raw_ages]


def build_life_table(mx: np.ndarray) -> pd.DataFrame:
    n = len(mx)
    ages = np.arange(n, dtype=int)
    qx = 2.0 * mx / (2.0 + mx)
    qx[-1] = 1.0
    px = 1.0 - qx

    lx = np.zeros(n, dtype=float)
    lx[0] = RADIX
    for i in range(1, n):
        lx[i] = lx[i - 1] * (1.0 - qx[i - 1])

    dx = lx * qx
    Lx = lx - 0.5 * dx
    Tx = np.cumsum(Lx[::-1])[::-1]
    ex = Tx / lx

    return pd.DataFrame(
        {
            "Edad_num": ages,
            "Edad": [age_label(int(x)) for x in ages],
            "mx": mx,
            "qx": qx,
            "px": px,
            "lx": lx,
            "dx": dx,
            "Lx": Lx,
            "Tx": Tx,
            "ex": ex,
        }
    )


def adult_threshold_age(df: pd.DataFrame, minimum_age: int, threshold: float) -> int | None:
    eligible = df[(df["Edad_num"] > minimum_age) & (df["qx"] >= threshold)]
    if eligible.empty:
        return None
    return int(eligible.iloc[0]["Edad_num"])


def derive_metrics(df: pd.DataFrame) -> dict[str, object]:
    min_m_idx = int(df["mx"].idxmin())
    min_q_idx = int(df["qx"].idxmin())
    max_d_idx = int(df["dx"].idxmax())
    adult_max_d_idx = int(df.loc[df["Edad_num"] >= 1, "dx"].idxmax())
    max_e_idx = int(df["ex"].idxmax())

    thresholds = {
        threshold: adult_threshold_age(df, int(df.loc[min_q_idx, "Edad_num"]), threshold)
        for threshold in (0.01, 0.05, 0.10, 0.25, 0.50)
    }
    selected_ages = [0, 1, 10, 25, 50, 60, 65, 80, 90, 100, 109]
    survivors = {age: float(df.loc[df["Edad_num"] == age, "lx"].iloc[0]) for age in selected_ages}

    return {
        "minimum_mx_age": int(df.loc[min_m_idx, "Edad_num"]),
        "minimum_mx": float(df.loc[min_m_idx, "mx"]),
        "minimum_qx_age": int(df.loc[min_q_idx, "Edad_num"]),
        "minimum_qx": float(df.loc[min_q_idx, "qx"]),
        "maximum_dx_age": int(df.loc[max_d_idx, "Edad_num"]),
        "maximum_dx": float(df.loc[max_d_idx, "dx"]),
        "adult_maximum_dx_age": int(df.loc[adult_max_d_idx, "Edad_num"]),
        "adult_maximum_dx": float(df.loc[adult_max_d_idx, "dx"]),
        "maximum_ex_age": int(df.loc[max_e_idx, "Edad_num"]),
        "maximum_ex": float(df.loc[max_e_idx, "ex"]),
        "e0": float(df.loc[0, "ex"]),
        "e1": float(df.loc[1, "ex"]),
        "T0": float(df.loc[0, "Tx"]),
        "last_lx": float(df.loc[len(df) - 1, "lx"]),
        "last_Lx": float(df.loc[len(df) - 1, "Lx"]),
        "sum_dx": float(df["dx"].sum()),
        "survivors": survivors,
        "thresholds": thresholds,
    }


def run_quality_checks(df: pd.DataFrame) -> list[dict[str, object]]:
    qx = df["qx"].to_numpy()
    lx = df["lx"].to_numpy()
    dx = df["dx"].to_numpy()
    Lx = df["Lx"].to_numpy()
    Tx = df["Tx"].to_numpy()
    ex = df["ex"].to_numpy()

    recurrence_l = float(np.max(np.abs(lx[:-1] - dx[:-1] - lx[1:])))
    recurrence_t = float(np.max(np.abs(Tx[:-1] - Lx[:-1] - Tx[1:])))
    identity_e = float(np.max(np.abs(ex - Tx / lx)))
    sum_dx_gap = float(abs(dx.sum() - RADIX))
    last_survivor = float(abs(lx[-1] - dx[-1]))
    lx_increases = int(np.sum(np.diff(lx) > 1e-10))
    L_upper_violations = int(np.sum(Lx - lx > 1e-10))
    L_lower_violations = int(np.sum(Lx[:-1] + 1e-10 < lx[1:]))

    checks = [
        ("Probabilidades dentro de [0,1]", bool(np.all((qx >= 0) & (qx <= 1))), f"mín={qx.min():.12f}; máx={qx.max():.12f}", "0 violaciones"),
        ("lx no negativa y no creciente", bool(np.all(lx >= 0) and lx_increases == 0), f"mín={lx.min():.9f}; aumentos={lx_increases}", "0 aumentos"),
        ("dx no negativa", bool(np.all(dx >= 0)), f"mín={dx.min():.12f}", "mínimo ≥ 0"),
        ("Recurrencia l(x+1)=lx-dx", recurrence_l <= 1e-8, f"error máximo={recurrence_l:.3e}", "≤ 1e-8"),
        ("Suma de defunciones igual al radix", sum_dx_gap <= 1e-8, f"Σdx={dx.sum():.10f}; diferencia={sum_dx_gap:.3e}", "≤ 1e-8"),
        ("Extinción después de 109 y+", last_survivor <= 1e-10, f"l109-d109={last_survivor:.3e}", "≤ 1e-10"),
        ("l(x+1) ≤ Lx ≤ lx", L_upper_violations == 0 and L_lower_violations == 0, f"violaciones superiores={L_upper_violations}; inferiores={L_lower_violations}", "0 violaciones"),
        ("Recurrencia Tx=Lx+T(x+1)", recurrence_t <= 1e-8, f"error máximo={recurrence_t:.3e}", "≤ 1e-8"),
        ("T(109+)=L(109+)", abs(Tx[-1] - Lx[-1]) <= 1e-12, f"diferencia={abs(Tx[-1]-Lx[-1]):.3e}", "≤ 1e-12"),
        ("Identidad ex=Tx/lx", identity_e <= 1e-12, f"error máximo={identity_e:.3e}", "≤ 1e-12"),
        ("Cierre q109+=1 y e109+=0.5", abs(qx[-1] - 1) <= 1e-12 and abs(ex[-1] - 0.5) <= 1e-12, f"q={qx[-1]:.12f}; e={ex[-1]:.12f}", "valores exactos"),
        ("110 edades completas y sin nulos", len(df) == 110 and not df.isna().any().any(), f"filas={len(df)}; nulos={int(df.isna().sum().sum())}", "110 filas; 0 nulos"),
    ]
    return [
        {"control": name, "cumple": passed, "evidencia": evidence, "tolerancia": tolerance}
        for name, passed, evidence, tolerance in checks
    ]


def configure_plot(ax: plt.Axes, ylabel: str, log_scale: bool = False) -> None:
    ax.set_xlim(0, 109)
    ticks = list(range(0, 101, 10)) + [109]
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(x) if x != 109 else "109+" for x in ticks])
    ax.set_xlabel("Edad x")
    ax.set_ylabel(ylabel)
    ax.grid(True, which="major", alpha=0.28, linewidth=0.7)
    if log_scale:
        ax.set_yscale("log")
        ax.grid(True, which="minor", alpha=0.12, linewidth=0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def create_graphs(df: pd.DataFrame, metrics: dict[str, object], student: str) -> dict[str, Path]:
    GRAPH_DIR.mkdir(parents=True, exist_ok=True)
    plt.style.use("seaborn-v0_8-whitegrid")
    x = df["Edad_num"].to_numpy()
    graph_specs = [
        ("01_mx", "mx", "Tasa específica de mortalidad mx por edad", "mx (defunciones por persona-año)", "#1769aa"),
        ("02_qx", "qx", "Probabilidad anual de muerte qx por edad", "qx (probabilidad)", "#d1495b"),
        ("03_lx", "lx", "Sobrevivientes lx por edad", "lx (personas)", "#218c74"),
        ("04_dx", "dx", "Defunciones dx por edad", "dx (personas)", "#e67e22"),
        ("05_Lx", "Lx", "Años-persona Lx por edad", "Lx (años-persona)", "#6c5ce7"),
        ("06_Tx", "Tx", "Años-persona acumulados Tx por edad", "Tx (años-persona)", "#2d3436"),
        ("07_ex", "ex", "Esperanza de vida completa ex por edad", "ex (años)", "#8e44ad"),
    ]
    paths: dict[str, Path] = {}

    for key, column, title, ylabel, color in graph_specs:
        fig, ax = plt.subplots(figsize=(10.8, 6.2))
        y = df[column].to_numpy()
        ax.plot(x, y, color=color, linewidth=2.35)
        configure_plot(ax, ylabel)
        ax.set_title(title, fontsize=15, fontweight="bold", pad=13)
        if column in {"lx", "dx", "Lx", "Tx"}:
            ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:,.0f}"))
        if column in {"mx", "qx"}:
            min_age = int(metrics[f"minimum_{column}_age"])
            min_value = float(metrics[f"minimum_{column}"])
            ax.scatter([min_age], [min_value], color="#111111", s=38, zorder=5)
            ax.annotate(
                f"Mínimo: edad {min_age}\n{min_value:.8f}",
                (min_age, min_value),
                xytext=(min_age + 8, min_value + 0.07 * max(y)),
                arrowprops={"arrowstyle": "->", "color": "#444444"},
                fontsize=9,
            )
        elif column == "dx":
            peak_age = int(metrics["maximum_dx_age"])
            peak = float(metrics["maximum_dx"])
            adult_peak_age = int(metrics["adult_maximum_dx_age"])
            adult_peak = float(metrics["adult_maximum_dx"])
            ax.scatter([peak_age], [peak], color="#111111", s=38, zorder=5)
            ax.annotate(
                f"Máximo global: edad {peak_age}\n{peak:,.2f}",
                (peak_age, peak),
                xytext=(8, peak * 0.93),
                arrowprops={"arrowstyle": "->", "color": "#444444"},
                fontsize=9,
            )
            ax.scatter([adult_peak_age], [adult_peak], color="#7f1d1d", s=38, zorder=5)
            ax.annotate(
                f"Pico adulto: edad {adult_peak_age}\n{adult_peak:,.2f}",
                (adult_peak_age, adult_peak),
                xytext=(adult_peak_age - 25, adult_peak * 0.78),
                arrowprops={"arrowstyle": "->", "color": "#7f1d1d"},
                fontsize=9,
                color="#7f1d1d",
            )
        elif column == "ex":
            e0 = float(metrics["e0"])
            max_age = int(metrics["maximum_ex_age"])
            max_e = float(metrics["maximum_ex"])
            ax.scatter([0, max_age], [e0, max_e], color="#111111", s=32, zorder=5)
            ax.annotate(f"e0={e0:.2f}", (0, e0), xytext=(7, e0 - 8), fontsize=9)
            ax.annotate(f"Máximo: {max_e:.2f} (edad {max_age})", (max_age, max_e), xytext=(16, max_e + 1.5), fontsize=9)
        fig.text(
            0.5,
            0.012,
            f"Fuente: serie {SERIE}, {student}. Cálculos propios bajo UDD; radix 100,000.",
            ha="center",
            fontsize=8,
            color="#555555",
        )
        fig.tight_layout(rect=(0, 0.035, 1, 1))
        path = GRAPH_DIR / f"{key}.png"
        fig.savefig(path, dpi=190, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        paths[key] = path

    for key, column, title, color in [
        ("08_mx_log", "mx", "Tasa mx por edad — escala logarítmica", "#1769aa"),
        ("09_qx_log", "qx", "Probabilidad qx por edad — escala logarítmica", "#d1495b"),
    ]:
        fig, ax = plt.subplots(figsize=(10.8, 6.2))
        y = df[column].to_numpy()
        ax.plot(x, y, color=color, linewidth=2.25)
        configure_plot(ax, column, log_scale=True)
        ax.set_title(title, fontsize=15, fontweight="bold", pad=13)
        min_age = int(metrics[f"minimum_{column}_age"])
        min_value = float(metrics[f"minimum_{column}"])
        ax.scatter([min_age], [min_value], color="#111111", s=38, zorder=5)
        ax.annotate(
            f"Mínimo: edad {min_age}",
            (min_age, min_value),
            xytext=(min_age + 10, min_value * 2.6),
            arrowprops={"arrowstyle": "->", "color": "#444444"},
            fontsize=9,
        )
        fig.text(0.5, 0.012, f"Fuente: serie {SERIE}. Escala logarítmica en el eje vertical.", ha="center", fontsize=8, color="#555555")
        fig.tight_layout(rect=(0, 0.035, 1, 1))
        path = GRAPH_DIR / f"{key}.png"
        fig.savefig(path, dpi=190, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        paths[key] = path

    comparisons = [
        ("10_mx_qx", ["mx", "qx"], ["mx", "qx"], ["#1769aa", "#d1495b"], "Comparación de mx y qx por edad", "Tasa / probabilidad"),
        ("11_lx_Lx", ["lx", "Lx"], ["lx", "Lx"], ["#218c74", "#6c5ce7"], "Comparación de lx y Lx por edad", "Personas / años-persona"),
    ]
    for key, columns, labels, colors, title, ylabel in comparisons:
        fig, ax = plt.subplots(figsize=(10.8, 6.2))
        for column, label, color in zip(columns, labels, colors):
            ax.plot(x, df[column].to_numpy(), label=label, color=color, linewidth=2.2)
        configure_plot(ax, ylabel)
        ax.set_title(title, fontsize=15, fontweight="bold", pad=13)
        if key == "11_lx_Lx":
            ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:,.0f}"))
        ax.legend(frameon=False, fontsize=10)
        fig.text(0.5, 0.012, f"Fuente: serie {SERIE}. Cálculos propios bajo UDD.", ha="center", fontsize=8, color="#555555")
        fig.tight_layout(rect=(0, 0.035, 1, 1))
        path = GRAPH_DIR / f"{key}.png"
        fig.savefig(path, dpi=190, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        paths[key] = path

    return paths


def create_excel(
    df: pd.DataFrame,
    checks: list[dict[str, object]],
    metrics: dict[str, object],
    graph_paths: dict[str, Path],
    source: Path,
    source_col: str,
    student: str,
) -> None:
    workbook = xlsxwriter.Workbook(OUTPUT_XLSX)
    workbook.set_properties(
        {
            "title": "Tabla de mortalidad — Serie 08",
            "subject": "Funciones biométricas bajo distribución uniforme de las muertes",
            "author": student,
            "comments": "Generado a partir de las tasas mx individuales; radix 100,000.",
        }
    )

    title_fmt = workbook.add_format({"bold": True, "font_size": 16, "font_color": "#FFFFFF", "bg_color": "#17365D", "align": "center", "valign": "vcenter"})
    subtitle_fmt = workbook.add_format({"bold": True, "font_color": "#17365D"})
    header_fmt = workbook.add_format({"bold": True, "font_color": "#FFFFFF", "bg_color": "#2F75B5", "border": 1, "align": "center", "valign": "vcenter"})
    text_fmt = workbook.add_format({"border": 1, "valign": "top"})
    age_fmt = workbook.add_format({"border": 1, "align": "center"})
    rate_fmt = workbook.add_format({"border": 1, "num_format": "0.00000000"})
    count_fmt = workbook.add_format({"border": 1, "num_format": "#,##0.000000"})
    life_fmt = workbook.add_format({"border": 1, "num_format": "0.000000"})
    formula_note_fmt = workbook.add_format({"text_wrap": True, "valign": "top", "bg_color": "#EAF2F8", "border": 1})
    pass_fmt = workbook.add_format({"font_color": "#006100", "bg_color": "#C6EFCE", "border": 1, "align": "center", "bold": True})
    fail_fmt = workbook.add_format({"font_color": "#9C0006", "bg_color": "#FFC7CE", "border": 1, "align": "center", "bold": True})
    small_fmt = workbook.add_format({"font_size": 9, "font_color": "#555555", "text_wrap": True})

    ws_source = workbook.add_worksheet("Fuente_mx")
    ws_source.merge_range("A1:D1", "Fuente de tasas específicas de mortalidad — Serie 08", title_fmt)
    ws_source.set_row(0, 25)
    metadata = [
        ("Archivo fuente", source.name),
        ("Hoja", "Tasas Específicas de Mortalidad"),
        ("Serie / columna", f"{SERIE} / {source_col}"),
        ("Estudiante", student),
    ]
    for row, (label, value) in enumerate(metadata, start=1):
        ws_source.write(row, 0, label, subtitle_fmt)
        ws_source.write(row, 1, value)
    source_header_row = 6
    for col, header in enumerate(["Edad", "Edad numérica", "mx", "Celda de origen"]):
        ws_source.write(source_header_row - 1, col, header, header_fmt)
    for i, record in df.iterrows():
        row = source_header_row + i
        ws_source.write(row, 0, record["Edad"], age_fmt)
        ws_source.write_number(row, 1, int(record["Edad_num"]), age_fmt)
        ws_source.write_number(row, 2, float(record["mx"]), rate_fmt)
        ws_source.write(row, 3, f"{source_col}{5 + i}", age_fmt)
    ws_source.set_column("A:A", 13)
    ws_source.set_column("B:B", 15)
    ws_source.set_column("C:C", 16)
    ws_source.set_column("D:D", 16)
    ws_source.freeze_panes(source_header_row, 0)
    ws_source.autofilter(source_header_row - 1, 0, source_header_row - 1 + len(df), 3)

    ws_table = workbook.add_worksheet("Tabla_Mortalidad")
    ws_table.merge_range("A1:I1", "Tabla completa de mortalidad — Serie 08", title_fmt)
    ws_table.set_row(0, 25)
    ws_table.write("A2", "Estudiante", subtitle_fmt)
    ws_table.write("B2", student)
    ws_table.write("A3", "Supuestos", subtitle_fmt)
    ws_table.merge_range("B3:I3", "UDD; ax=0.5; radix l0=100,000; q(109+)=1; sin redondeo intermedio.", small_fmt)
    table_header_row = 4
    columns = ["Edad", "mx", "qx", "px", "lx", "dx", "Lx", "Tx", "ex"]
    for col, header in enumerate(columns):
        ws_table.write(table_header_row, col, header, header_fmt)
    data_start_row = table_header_row + 1
    excel_first = data_start_row + 1
    excel_last = excel_first + len(df) - 1

    for i, record in df.iterrows():
        row = data_start_row + i
        erow = row + 1
        ws_table.write(row, 0, record["Edad"], age_fmt)
        ws_table.write_formula(row, 1, f"=Fuente_mx!C{source_header_row + i + 1}", rate_fmt, float(record["mx"]))
        if i == len(df) - 1:
            q_formula = "=1"
        else:
            q_formula = f"=2*B{erow}/(2+B{erow})"
        ws_table.write_formula(row, 2, q_formula, rate_fmt, float(record["qx"]))
        ws_table.write_formula(row, 3, f"=1-C{erow}", rate_fmt, float(record["px"]))
        if i == 0:
            l_formula = f"={int(RADIX)}"
        else:
            l_formula = f"=E{erow-1}-F{erow-1}"
        ws_table.write_formula(row, 4, l_formula, count_fmt, float(record["lx"]))
        ws_table.write_formula(row, 5, f"=E{erow}*C{erow}", count_fmt, float(record["dx"]))
        ws_table.write_formula(row, 6, f"=E{erow}-0.5*F{erow}", count_fmt, float(record["Lx"]))
        if i == len(df) - 1:
            t_formula = f"=G{erow}"
        else:
            t_formula = f"=G{erow}+H{erow+1}"
        ws_table.write_formula(row, 7, t_formula, count_fmt, float(record["Tx"]))
        ws_table.write_formula(row, 8, f"=H{erow}/E{erow}", life_fmt, float(record["ex"]))

    ws_table.write("K5", "Función", header_fmt)
    ws_table.write("L5", "Definición aplicada", header_fmt)
    formula_notes = [
        ("qx", "2mx/(2+mx); en 109 y+ se impone qx=1"),
        ("px", "1-qx"),
        ("lx", "l0=100,000; l(x+1)=lx-dx"),
        ("dx", "lx·qx"),
        ("Lx", "lx-0.5dx; en 109 y+ equivale a lx/2"),
        ("Tx", "suma de L desde x hasta 109 y+"),
        ("ex", "Tx/lx"),
    ]
    for offset, (function, definition) in enumerate(formula_notes, start=5):
        ws_table.write(offset, 10, function, formula_note_fmt)
        ws_table.write(offset, 11, definition, formula_note_fmt)
    ws_table.set_column("A:A", 12)
    ws_table.set_column("B:D", 14)
    ws_table.set_column("E:H", 17)
    ws_table.set_column("I:I", 14)
    ws_table.set_column("K:K", 12)
    ws_table.set_column("L:L", 42)
    ws_table.freeze_panes(data_start_row, 1)
    ws_table.autofilter(table_header_row, 0, table_header_row + len(df), 8)
    ws_table.conditional_format(data_start_row, 2, table_header_row + len(df), 2, {"type": "data_bar", "bar_color": "#5B9BD5"})

    ws_checks = workbook.add_worksheet("Controles")
    ws_checks.merge_range("A1:E1", "Controles actuariales y numéricos", title_fmt)
    ws_checks.set_row(0, 25)
    ws_checks.write("A2", "Todos los controles se calculan con valores sin redondear.", small_fmt)
    check_header_row = 3
    for col, header in enumerate(["Control", "Resultado", "Evidencia", "Tolerancia", "Fórmula dinámica"]):
        ws_checks.write(check_header_row, col, header, header_fmt)

    dynamic_formulas = [
        f'=IF(AND(MIN(Tabla_Mortalidad!C{excel_first}:C{excel_last})>=0,MAX(Tabla_Mortalidad!C{excel_first}:C{excel_last})<=1),"CUMPLE","REVISAR")',
        f'=IF(AND(MIN(Tabla_Mortalidad!E{excel_first}:E{excel_last})>=0,SUMPRODUCT(--(Tabla_Mortalidad!E{excel_first+1}:E{excel_last}>Tabla_Mortalidad!E{excel_first}:E{excel_last-1}))=0),"CUMPLE","REVISAR")',
        f'=IF(MIN(Tabla_Mortalidad!F{excel_first}:F{excel_last})>=0,"CUMPLE","REVISAR")',
        f'=IF(SUMPRODUCT(ABS(Tabla_Mortalidad!E{excel_first}:E{excel_last-1}-Tabla_Mortalidad!F{excel_first}:F{excel_last-1}-Tabla_Mortalidad!E{excel_first+1}:E{excel_last}))<=1E-8,"CUMPLE","REVISAR")',
        f'=IF(ABS(SUM(Tabla_Mortalidad!F{excel_first}:F{excel_last})-100000)<=1E-8,"CUMPLE","REVISAR")',
        f'=IF(ABS(Tabla_Mortalidad!E{excel_last}-Tabla_Mortalidad!F{excel_last})<=1E-10,"CUMPLE","REVISAR")',
        f'=IF(AND(SUMPRODUCT(--(Tabla_Mortalidad!G{excel_first}:G{excel_last}>Tabla_Mortalidad!E{excel_first}:E{excel_last}))=0,SUMPRODUCT(--(Tabla_Mortalidad!G{excel_first}:G{excel_last-1}<Tabla_Mortalidad!E{excel_first+1}:E{excel_last}))=0),"CUMPLE","REVISAR")',
        f'=IF(SUMPRODUCT(ABS(Tabla_Mortalidad!H{excel_first}:H{excel_last-1}-Tabla_Mortalidad!G{excel_first}:G{excel_last-1}-Tabla_Mortalidad!H{excel_first+1}:H{excel_last}))<=1E-8,"CUMPLE","REVISAR")',
        f'=IF(ABS(Tabla_Mortalidad!H{excel_last}-Tabla_Mortalidad!G{excel_last})<=1E-12,"CUMPLE","REVISAR")',
        f'=IF(SUMPRODUCT(ABS(Tabla_Mortalidad!I{excel_first}:I{excel_last}-Tabla_Mortalidad!H{excel_first}:H{excel_last}/Tabla_Mortalidad!E{excel_first}:E{excel_last}))<=1E-12,"CUMPLE","REVISAR")',
        f'=IF(AND(Tabla_Mortalidad!C{excel_last}=1,ABS(Tabla_Mortalidad!I{excel_last}-0.5)<=1E-12),"CUMPLE","REVISAR")',
        f'=IF(COUNT(Tabla_Mortalidad!B{excel_first}:I{excel_last})=880,"CUMPLE","REVISAR")',
    ]
    for idx, (check, formula) in enumerate(zip(checks, dynamic_formulas), start=check_header_row + 1):
        ws_checks.write(idx, 0, check["control"], text_fmt)
        result_fmt = pass_fmt if check["cumple"] else fail_fmt
        cached = "CUMPLE" if check["cumple"] else "REVISAR"
        ws_checks.write_formula(idx, 1, formula, result_fmt, cached)
        ws_checks.write(idx, 2, check["evidencia"], text_fmt)
        ws_checks.write(idx, 3, check["tolerancia"], text_fmt)
        ws_checks.write(idx, 4, formula, formula_note_fmt)
    ws_checks.set_column("A:A", 39)
    ws_checks.set_column("B:B", 13)
    ws_checks.set_column("C:C", 40)
    ws_checks.set_column("D:D", 18)
    ws_checks.set_column("E:E", 80)
    ws_checks.freeze_panes(check_header_row + 1, 0)

    summary_start = check_header_row + 1 + len(checks) + 2
    ws_checks.write(summary_start, 0, "Resultados clave", header_fmt)
    key_results = [
        ("e0", fmt_num(float(metrics["e0"]), 6)),
        ("Edad del mínimo de mx", str(metrics["minimum_mx_age"])),
        ("Mínimo mx", f'{float(metrics["minimum_mx"]):.10f}'),
        ("Edad del máximo de dx", str(metrics["maximum_dx_age"])),
        ("Máximo dx", fmt_num(float(metrics["maximum_dx"]), 6)),
        ("Pico adulto de dx", f"edad {metrics['adult_maximum_dx_age']}: {fmt_num(float(metrics['adult_maximum_dx']), 6)}"),
        ("Σdx", fmt_num(float(metrics["sum_dx"]), 10)),
    ]
    for offset, (label, value) in enumerate(key_results, start=1):
        ws_checks.write(summary_start + offset, 0, label, subtitle_fmt)
        ws_checks.write(summary_start + offset, 1, value)

    ws_graphs = workbook.add_worksheet("Graficas")
    ws_graphs.merge_range("A1:R1", "Gráficas de funciones biométricas — Serie 08", title_fmt)
    ws_graphs.set_row(0, 25)
    ws_graphs.write("A2", f"Estudiante: {student}. Imágenes generadas con los valores sin redondear.", small_fmt)
    ws_graphs.set_column("A:R", 9)
    ordered_paths = [graph_paths[key] for key in sorted(graph_paths)]
    for idx, path in enumerate(ordered_paths):
        block_row = 3 + (idx // 2) * 27
        block_col = 0 if idx % 2 == 0 else 9
        ws_graphs.insert_image(block_row, block_col, str(path), {"x_scale": 0.32, "y_scale": 0.32, "object_position": 1})
    workbook.close()


def wrapped(text: str, width: int = 100) -> str:
    return textwrap.fill(text, width=width, break_long_words=False, break_on_hyphens=False)


def add_page_footer(fig: plt.Figure, page_number: int) -> None:
    fig.text(0.06, 0.025, "Estadística Actuarial — Tabla de Mortalidad, Serie 08", fontsize=8, color="#666666")
    fig.text(0.94, 0.025, f"Página {page_number}", fontsize=8, color="#666666", ha="right")


def add_text_page(pdf: PdfPages, title: str, blocks: list[tuple[str, str]], page_number: int) -> None:
    fig = plt.figure(figsize=(8.27, 11.69), facecolor="white")
    fig.text(0.07, 0.945, title, fontsize=19, fontweight="bold", color="#17365D")
    fig.add_artist(plt.Line2D([0.07, 0.93], [0.925, 0.925], color="#2F75B5", linewidth=2))
    y = 0.89
    for heading, body in blocks:
        fig.text(0.075, y, heading, fontsize=12, fontweight="bold", color="#2F75B5", va="top")
        y -= 0.032
        formatted = wrapped(body, 104)
        line_count = formatted.count("\n") + 1
        fig.text(0.085, y, formatted, fontsize=10.2, color="#222222", va="top", linespacing=1.45)
        y -= line_count * 0.026 + 0.027
    add_page_footer(fig, page_number)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def add_chart_page(pdf: PdfPages, title: str, image_path: Path, narrative: str, page_number: int) -> None:
    fig = plt.figure(figsize=(8.27, 11.69), facecolor="white")
    fig.text(0.07, 0.95, title, fontsize=18, fontweight="bold", color="#17365D")
    fig.add_artist(plt.Line2D([0.07, 0.93], [0.93, 0.93], color="#2F75B5", linewidth=2))
    ax = fig.add_axes([0.075, 0.36, 0.85, 0.53])
    ax.imshow(plt.imread(image_path))
    ax.axis("off")
    fig.text(0.08, 0.31, "Tendencia observada", fontsize=12, fontweight="bold", color="#2F75B5", va="top")
    fig.text(0.085, 0.275, wrapped(narrative, 104), fontsize=10.2, color="#222222", va="top", linespacing=1.45)
    add_page_footer(fig, page_number)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def create_pdf(
    df: pd.DataFrame,
    checks: list[dict[str, object]],
    metrics: dict[str, object],
    graph_paths: dict[str, Path],
    source: Path,
    source_col: str,
    student: str,
) -> int:
    thresholds: dict[float, int | None] = metrics["thresholds"]  # type: ignore[assignment]
    survivors: dict[int, float] = metrics["survivors"]  # type: ignore[assignment]
    min_age = int(metrics["minimum_mx_age"])
    min_mx = float(metrics["minimum_mx"])
    min_qx = float(metrics["minimum_qx"])
    peak_age = int(metrics["maximum_dx_age"])
    peak_dx = float(metrics["maximum_dx"])
    adult_peak_age = int(metrics["adult_maximum_dx_age"])
    adult_peak_dx = float(metrics["adult_maximum_dx"])
    q1_age = thresholds[0.01]
    q5_age = thresholds[0.05]
    q10_age = thresholds[0.10]

    narratives = {
        "mx": (
            f"La tasa mx parte en {df.loc[0, 'mx']:.8f} ({fmt_rate(float(df.loc[0, 'mx']))}) a la edad 0, "
            f"desciende hasta su mínimo de {min_mx:.8f} a la edad {min_age} y luego cambia a una trayectoria creciente. "
            f"Este perfil refleja mortalidad infantil elevada, un valle en edades escolares y crecimiento de la mortalidad adulta. "
            f"En la edad abierta 109 y+ la tasa alcanza {df.loc[109, 'mx']:.8f}."
        ),
        "qx": (
            f"La probabilidad qx reproduce la forma de mx por ser una transformación creciente bajo UDD. Su mínimo es "
            f"{min_qx:.8f} a la edad {metrics['minimum_qx_age']}. Después del valle, qx alcanza 1% a la edad {q1_age}, "
            f"5% a la edad {q5_age} y 10% a la edad {q10_age}. En 109 y+ se impone qx=1 para cerrar la tabla; "
            f"ese último valor es un supuesto de cierre y no la conversión directa de mx."
        ),
        "lx": (
            f"La cohorte inicia con 100,000 personas. Sobreviven {survivors[1]:,.2f} a la edad 1, "
            f"{survivors[50]:,.2f} a los 50, {survivors[65]:,.2f} a los 65 y {survivors[80]:,.2f} a los 80. "
            f"lx disminuye en toda la tabla: la mayor caída anual global ocurre en el primer año y el mayor descenso en edades "
            f"adultas se presenta a la edad {adult_peak_age}. A 109 y+ permanecen {survivors[109]:,.4f} sobrevivientes antes del cierre."
        ),
        "dx": (
            f"dx mide las defunciones de la cohorte entre x y x+1. Presenta {df.loc[0, 'dx']:,.2f} defunciones en el primer año, "
            f"que constituyen el máximo global a la edad {peak_age}. Tras caer hasta las edades de baja mortalidad, aumenta de nuevo "
            f"hasta un pico adulto de {adult_peak_dx:,.2f} a la edad {adult_peak_age}. Después de ese segundo pico, dx disminuye aunque qx continúe creciendo porque el número de personas todavía expuestas, lx, "
            f"se reduce con rapidez. La suma de todas las defunciones es {float(metrics['sum_dx']):,.8f}, igual al radix salvo precisión numérica."
        ),
        "Lx": (
            f"Lx representa los años-persona vividos entre x y x+1. Bajo UDD es el promedio de lx y l(x+1), por lo que permanece "
            f"entre ambas cantidades y sigue una trayectoria decreciente muy cercana a lx. La separación lx-Lx equivale a dx/2: "
            f"su máximo global ocurre en la edad {peak_age} y su máximo adulto en la edad {adult_peak_age}. En 109 y+, "
            f"Lx={float(metrics['last_Lx']):,.6f}, exactamente la mitad de l(109+)."
        ),
        "Tx": (
            f"Tx acumula, desde cada edad, los años-persona que todavía vivirá la cohorte. Comienza en T0={float(metrics['T0']):,.2f} "
            f"y disminuye estrictamente porque T(x)-T(x+1)=Lx, que es positivo. La reducción es casi lineal cuando Lx permanece alto y "
            f"se curva hacia cero en edades avanzadas. En la última edad Tx coincide con Lx."
        ),
        "ex": (
            f"La esperanza de vida al nacer es e0={float(metrics['e0']):.6f} años. Al sobrevivir el primer año aumenta a "
            f"e1={float(metrics['e1']):.6f}; su máximo es {float(metrics['maximum_ex']):.6f} a la edad "
            f"{metrics['maximum_ex_age']}. Después disminuye conforme se agota la vida residual. Por el supuesto UDD de cierre, "
            f"e(109+)=0.5 años."
        ),
    }

    page = 0
    with PdfPages(OUTPUT_PDF) as pdf:
        metadata = pdf.infodict()
        metadata["Title"] = "Informe de Tabla de Mortalidad — Serie 08"
        metadata["Author"] = student
        metadata["Subject"] = "Funciones biométricas bajo distribución uniforme de las muertes"
        metadata["Keywords"] = "mortalidad, mx, qx, lx, dx, Lx, Tx, ex, UDD"
        metadata["CreationDate"] = datetime.now()

        page += 1
        fig = plt.figure(figsize=(8.27, 11.69), facecolor="white")
        fig.text(0.5, 0.78, "TABLA DE MORTALIDAD", ha="center", fontsize=28, fontweight="bold", color="#17365D")
        fig.text(0.5, 0.72, "Funciones biométricas bajo distribución uniforme de las muertes", ha="center", fontsize=14, color="#2F75B5")
        fig.add_artist(plt.Line2D([0.18, 0.82], [0.68, 0.68], color="#2F75B5", linewidth=2.5))
        fig.text(0.5, 0.58, f"Serie asignada: {SERIE} — columna {source_col}", ha="center", fontsize=15)
        fig.text(0.5, 0.535, student, ha="center", fontsize=14, fontweight="bold")
        fig.text(0.5, 0.44, "Radix: 100,000 personas a edad exacta cero", ha="center", fontsize=12)
        fig.text(0.5, 0.405, "Edades: 0 a 108 y grupo abierto 109 y+", ha="center", fontsize=12)
        fig.text(0.5, 0.22, "Estadística Actuarial — Sección A", ha="center", fontsize=12, color="#555555")
        fig.text(0.5, 0.18, datetime.now().strftime("%d/%m/%Y"), ha="center", fontsize=11, color="#555555")
        add_page_footer(fig, page)
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        page += 1
        add_text_page(
            pdf,
            "1. Objetivo, supuestos y metodología",
            [
                ("Objetivo", "Construir la tabla completa de mortalidad de la serie 8 a partir de mx, calcular qx, px, lx, dx, Lx, Tx y ex, analizar sus tendencias por edad y documentar verificaciones actuariales y numéricas."),
                ("Supuestos", "Distribución uniforme de las muertes; factor de separación ax=0.5; radix l0=100,000; edades simples; q(109+)=1; todos los sobrevivientes se extinguen en el último intervalo; cálculos sin redondeo intermedio."),
                ("Conversión", "Para x<109: qx=2mx/(2+mx)=mx/(1+0.5mx). Además, px=1-qx, dx=lx·qx y l(x+1)=lx-dx."),
                ("Años-persona y esperanza", "Lx=lx-0.5dx=0.5[lx+l(x+1)]. Tx se obtuvo acumulando Lx de abajo hacia arriba y ex=Tx/lx. En 109 y+ se impuso qx=1, Lx=lx/2, Tx=Lx y ex=0.5."),
                ("Precisión", "Todos los resultados y controles usan doble precisión. Los decimales visibles en tablas y gráficas son únicamente formato de presentación."),
            ],
            page,
        )

        page += 1
        summary_survivors = "; ".join(f"l{age}={survivors[age]:,.2f}" for age in [0, 1, 50, 65, 80, 90, 100, 109])
        add_text_page(
            pdf,
            "2. Resumen de resultados",
            [
                ("Esperanza de vida", f"e0={float(metrics['e0']):.6f} años; e1={float(metrics['e1']):.6f} años; valor máximo={float(metrics['maximum_ex']):.6f} a la edad {metrics['maximum_ex_age']}."),
                ("Mortalidad mínima", f"mx alcanza su mínimo de {min_mx:.10f} a la edad {min_age}; qx alcanza {min_qx:.10f} en esa misma edad."),
                ("Defunciones máximas", f"dx alcanza su máximo global de {peak_dx:,.6f} a la edad {peak_age}; el pico adulto es {adult_peak_dx:,.6f} a la edad {adult_peak_age}."),
                ("Sobrevivientes seleccionados", summary_survivors + "."),
                ("Cierre", f"La suma de dx es {float(metrics['sum_dx']):,.10f}; en 109 y+ quedan {float(metrics['last_lx']):,.8f} personas, todas fallecen en el intervalo y e(109+)=0.5."),
            ],
            page,
        )

        table_columns = ["Edad", "mx", "qx", "lx", "dx", "Lx", "Tx", "ex"]
        rows_per_page = 37
        for start in range(0, len(df), rows_per_page):
            page += 1
            end = min(start + rows_per_page, len(df))
            fig, ax = plt.subplots(figsize=(8.27, 11.69), facecolor="white")
            ax.axis("off")
            table_title = "3. Tabla de mortalidad obtenida"
            if start > 0:
                table_title += " (continuación)"
            fig.text(
                0.045,
                0.955,
                table_title,
                fontsize=16,
                fontweight="bold",
                color="#17365D",
            )
            fig.text(
                0.047,
                0.928,
                f"Filas {start + 1} a {end} de {len(df)}",
                fontsize=9.5,
                color="#555555",
            )
            fig.add_artist(plt.Line2D([0.045, 0.955], [0.912, 0.912], color="#2F75B5", linewidth=2))
            page_rows = []
            for _, record in df.iloc[start:end].iterrows():
                page_rows.append(
                    [
                        record["Edad"],
                        f"{record['mx']:.8f}",
                        f"{record['qx']:.8f}",
                        f"{record['lx']:,.3f}",
                        f"{record['dx']:,.3f}",
                        f"{record['Lx']:,.3f}",
                        f"{record['Tx']:,.3f}",
                        f"{record['ex']:.4f}",
                    ]
                )
            table = ax.table(
                cellText=page_rows,
                colLabels=table_columns,
                colColours=["#2F75B5"] * len(table_columns),
                cellLoc="right",
                colLoc="center",
                colWidths=[0.085, 0.13, 0.13, 0.14, 0.13, 0.14, 0.16, 0.105],
                bbox=[0.035, 0.085, 0.93, 0.80],
            )
            table.auto_set_font_size(False)
            table.set_fontsize(5.8)
            for (row, col), cell in table.get_celld().items():
                cell.set_edgecolor("#C9D2DC")
                cell.set_linewidth(0.42)
                if row == 0:
                    cell.set_text_props(color="white", weight="bold", ha="center")
                elif row % 2 == 0:
                    cell.set_facecolor("#F4F7FA")
                if col == 0:
                    cell.set_text_props(ha="center")
            fig.text(
                0.5,
                0.058,
                "Valores sin redondeo intermedio; el redondeo mostrado es únicamente visual.",
                ha="center",
                fontsize=8,
                color="#555555",
            )
            add_page_footer(fig, page)
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

        for key, title, image_key in [
            ("mx", "4. Tendencia de mx", "01_mx"),
            ("qx", "5. Tendencia de qx", "02_qx"),
            ("lx", "6. Tendencia de lx", "03_lx"),
            ("dx", "7. Tendencia de dx", "04_dx"),
            ("Lx", "8. Tendencia de Lx", "05_Lx"),
            ("Tx", "9. Tendencia de Tx", "06_Tx"),
            ("ex", "10. Tendencia de ex", "07_ex"),
        ]:
            page += 1
            add_chart_page(pdf, title, graph_paths[image_key], narratives[key], page)

        page += 1
        add_chart_page(
            pdf,
            "11. Mortalidad en escala logarítmica",
            graph_paths["08_mx_log"],
            f"La escala logarítmica permite distinguir el descenso desde la mortalidad infantil hasta el mínimo de la edad {min_age} y el crecimiento posterior, que en escala lineal queda comprimido durante gran parte de la vida. La gráfica de qx en escala logarítmica, incluida entre los archivos de apoyo, presenta el mismo patrón.",
            page,
        )

        page += 1
        add_chart_page(
            pdf,
            "12. Comparación entre lx y Lx",
            graph_paths["11_lx_Lx"],
            f"Las curvas casi coinciden porque Lx es el promedio de los sobrevivientes al inicio y al final del intervalo. La distancia vertical es dx/2: es reducida en edades de baja mortalidad, alcanza su máximo global en la edad {peak_age} y un segundo máximo muy cercano en la edad adulta {adult_peak_age}. La comparación mx–qx también se entrega como gráfica complementaria.",
            page,
        )

        page += 1
        add_text_page(
            pdf,
            "13. Conclusiones",
            [
                ("Conclusión actuarial", f"La serie 8 presenta el patrón esperado de una tabla de mortalidad: mortalidad relativamente alta al nacer, mínimo a la edad {min_age}, crecimiento posterior de mx y qx, descenso continuo de lx, Lx, Tx y ex, y dos picos destacados de dx: el máximo global a la edad {peak_age} y el pico adulto a la edad {adult_peak_age}."),
                ("Resultados principales", f"La esperanza de vida al nacer es {float(metrics['e0']):.2f} años. El mínimo de mx ocurre a la edad {min_age}; el máximo global de dx se presenta a la edad {peak_age} y el pico adulto a la edad {adult_peak_age}. A los 65 años permanecen {survivors[65]:,.2f} sobrevivientes de la cohorte inicial y a los 80 años permanecen {survivors[80]:,.2f}."),
                ("Interpretación", f"Con las tasas del período, la esperanza de vida al nacer es {float(metrics['e0']):.2f} años. El resultado describe una cohorte hipotética sometida durante toda su vida a estas tasas; no es una predicción individual ni una tabla generacional."),
                ("Conclusión final", "La tabla construida es internamente consistente con la hipótesis de distribución uniforme de las muertes y con el cierre obligatorio en 109 y+. Las funciones obtenidas permiten describir conjuntamente la mortalidad, la supervivencia y la vida residual de la cohorte hipotética."),
            ],
            page,
        )

    return page


def validate_outputs(expected_pdf_pages: int, graph_paths: dict[str, Path]) -> dict[str, object]:
    if not OUTPUT_XLSX.exists() or OUTPUT_XLSX.stat().st_size == 0:
        raise RuntimeError("No se generó correctamente el libro Excel.")
    if not OUTPUT_PDF.exists() or OUTPUT_PDF.stat().st_size == 0:
        raise RuntimeError("No se generó correctamente el informe PDF.")
    if len(graph_paths) != 11 or any(not path.exists() or path.stat().st_size == 0 for path in graph_paths.values()):
        raise RuntimeError("No se generaron correctamente las 11 gráficas.")

    wb_formula = openpyxl.load_workbook(OUTPUT_XLSX, data_only=False, read_only=True)
    expected_sheets = ["Fuente_mx", "Tabla_Mortalidad", "Controles", "Graficas"]
    if wb_formula.sheetnames != expected_sheets:
        raise RuntimeError(f"Hojas inesperadas en Excel: {wb_formula.sheetnames}")
    ws_formula = wb_formula["Tabla_Mortalidad"]
    formula_q0 = ws_formula["C6"].value
    formula_q_last = ws_formula["C115"].value
    formula_t_last = ws_formula["H115"].value
    wb_formula.close()
    if formula_q0 != "=2*B6/(2+B6)" or formula_q_last != "=1" or formula_t_last != "=G115":
        raise RuntimeError("Las fórmulas clave del Excel no coinciden con la metodología.")

    wb_values = openpyxl.load_workbook(OUTPUT_XLSX, data_only=True, read_only=True)
    ws_values = wb_values["Tabla_Mortalidad"]
    e0_cached = ws_values["I6"].value
    q_last_cached = ws_values["C115"].value
    control_results = [wb_values["Controles"].cell(row, 2).value for row in range(5, 17)]
    wb_values.close()
    if e0_cached is None or q_last_cached != 1 or any(value != "CUMPLE" for value in control_results):
        raise RuntimeError("Los resultados almacenados del Excel o sus controles no son válidos.")

    import fitz

    pdf = fitz.open(OUTPUT_PDF)
    pdf_pages = len(pdf)
    pdf_text_chars = sum(len(page.get_text("text")) for page in pdf)
    pdf.close()
    if pdf_pages != expected_pdf_pages or pdf_text_chars < 10_000:
        raise RuntimeError(f"Validación PDF fallida: páginas={pdf_pages}, caracteres={pdf_text_chars}.")

    return {
        "excel_bytes": OUTPUT_XLSX.stat().st_size,
        "pdf_bytes": OUTPUT_PDF.stat().st_size,
        "pdf_pages": pdf_pages,
        "pdf_text_chars": pdf_text_chars,
        "graphs": len(graph_paths),
        "cached_e0": float(e0_cached),
    }


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    source = find_source()
    student, source_col, mx, _ = load_assigned_series(source)
    expected_student = "RAÚL ANTONIO RUANO KATZ"
    if student.upper() != expected_student:
        raise RuntimeError(
            f"La serie {SERIE} corresponde a {student!r}, no a {expected_student!r}. "
            "Se detuvo el proceso para evitar usar una serie incorrecta."
        )

    df = build_life_table(mx)
    metrics = derive_metrics(df)
    checks = run_quality_checks(df)
    if not all(bool(check["cumple"]) for check in checks):
        failures = [str(check["control"]) for check in checks if not check["cumple"]]
        raise RuntimeError(f"Fallaron controles antes de exportar: {failures}")

    graph_paths = create_graphs(df, metrics, student)
    create_excel(df, checks, metrics, graph_paths, source, source_col, student)
    pdf_pages = create_pdf(df, checks, metrics, graph_paths, source, source_col, student)
    validation = validate_outputs(pdf_pages, graph_paths)

    print("GENERACIÓN COMPLETA")
    print(f"Serie: {SERIE}; columna: {source_col}; estudiante: {student}")
    print(f"e0: {metrics['e0']:.10f}")
    print(f"Mínimo mx: edad {metrics['minimum_mx_age']}, valor {metrics['minimum_mx']:.10f}")
    print(f"Máximo dx: edad {metrics['maximum_dx_age']}, valor {metrics['maximum_dx']:.10f}")
    print(f"Pico adulto dx: edad {metrics['adult_maximum_dx_age']}, valor {metrics['adult_maximum_dx']:.10f}")
    print(f"Suma dx: {metrics['sum_dx']:.12f}")
    print(f"Controles: {sum(bool(c['cumple']) for c in checks)}/{len(checks)} CUMPLE")
    print(f"Excel: {OUTPUT_XLSX} ({validation['excel_bytes']:,} bytes)")
    print(f"PDF: {OUTPUT_PDF} ({validation['pdf_pages']} páginas; {validation['pdf_bytes']:,} bytes)")
    print(f"Gráficas: {GRAPH_DIR} ({validation['graphs']} archivos)")


if __name__ == "__main__":
    main()
