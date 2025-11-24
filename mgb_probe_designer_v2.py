#!/usr/bin/env python3
"""
MGB Probe Multiplex Designer v2.0
Academic Release — Single file, fully self-contained
GitHub: https://github.com/yourusername/MGB-Probe-Designer
"""

import sys
import os
import re
import numpy as np
from dataclasses import dataclass
from typing import List, Dict, Optional
from enum import Enum
from itertools import combinations
from collections import defaultdict

# ------------------- Dependencies -------------------
try:
    from Bio import SeqIO
    import pandas as pd
    from PyQt6.QtWidgets import *
    from PyQt6.QtCore import Qt, QThread, pyqtSignal
    from PyQt6.QtGui import QFont, QPalette, QColor
except ImportError as e:
    print("Missing packages. Run:\npip install PyQt6 pandas openpyxl biopython numpy")
    sys.exit(1)


# ------------------- Core Classes -------------------
class Fluorophore(Enum):
    FAM = "FAM"
    VIC = "VIC"
    NED = "NED"
    ROX = "ROX"
    CY5 = "CY5"


@dataclass
class MGBProbe:
    sequence: str
    fluorophore: Optional[Fluorophore] = None
    target_allele: str = ""
    tm: float = 0.0
    dg: float = 0.0
    snp_position: int = 0
    gc_content: float = 0.0
    length: int = 0
    target_id: str = ""


class ThermodynamicsCalculator:
    NN_DH = {'AA': -7.9, 'TT': -7.9, 'AT': -7.2, 'TA': -7.2, 'CA': -8.5, 'TG': -8.5,
                      'GT': -8.4, 'AC': -8.4, 'CT': -7.8, 'AG': -7.8, 'GA': -8.2, 'TC': -8.2,
                      'CG': -10.6, 'GC': -9.8, 'GG': -8.0, 'CC': -8.0}
    NN_DS = {'AA': -22.2, 'TT': -22.2, 'AT': -20.4, 'TA': -21.3, 'CA': -22.7, 'TG': -22.7,
                      'GT': -22.4, 'AC': -22.4, 'CT': -21.0, 'AG': -21.0, 'GA': -22.2, 'TC': -22.2,
                      'CG': -27.2, 'GC': -24.4, 'GG': -19.9, 'CC': -19.9}
    MGB_BOOST = 15.0

    def calculate_tm(self, seq: str) -> float:
        seq = seq.upper()
        if len(seq) < 2: return 0.0
        dh = ds = 0.0
        for i in range(len(seq)-1):
            pair = seq[i:i+2]
            dh += self.NN_DH.get(pair, -8.0)
            ds += self.NN_DS.get(pair, -22.0)
        dh += 0.2; ds += -5.7
        ds += 0.368 * (len(seq)-1) * np.log(0.05)
        tm = (1000 * dh) / (ds + 1.987 * np.log(0.00025)) - 273.15
        return round(tm + self.MGB_BOOST, 2)

    def calculate_dg(self, seq: str, temp: float = 60.0) -> float:
        seq = seq.upper()
        dh = ds = 0.0
        for i in range(len(seq)-1):
            pair = seq[i:i+2]
            dh += self.NN_DH.get(pair, -8.0)
            ds += self.NN_DS.get(pair, -22.0)
        dg = dh - (temp + 273.15) * ds / 1000
        return round(dg, 2)

    def gc_content(self, seq: str) -> float:
        s = seq.upper()
        return round(100 * (s.count('G') + s.count('C')) / len(s), 1) if s else 0.0


class MGBProbeDesigner:
    def __init__(self, thermo: ThermodynamicsCalculator):
        self.thermo = thermo

    def design_for_snp(self, amplicon: str, snp_pos: int, ref: str, alt: str, target_id: str) -> List[MGBProbe]:
        seq = amplicon.upper()
        candidates = []
        for dist in range(1, 7):
            for length in range(15, 23):
                start = snp_pos - (length - dist)
                end = snp_pos + dist
                if start < 0 or end > len(seq): continue
                for allele, base in [("REF", ref), ("ALT", alt)]:
                    probe_seq = seq[start:snp_pos] + base + seq[snp_pos+1:end]
                    if re.search(r'AAAAA|TTTTT|CCCCC|GGGGG', probe_seq): continue
                    tm = self.thermo.calculate_tm(probe_seq)
                    if 64 <= tm <= 73:
                        candidates.append(MGBProbe(
                            sequence=probe_seq,
                            target_allele=allele,
                            tm=tm,
                            dg=self.thermo.calculate_dg(probe_seq),
                            snp_position=dist,
                            gc_content=self.thermo.gc_content(probe_seq),
                            length=length,
                            target_id=target_id
                        ))
        candidates.sort(key=lambda p: abs(p.tm - 68))
        best = []
        seen = set()
        for p in candidates:
            key = (p.target_id, p.target_allele)
            if key not in seen:
                best.append(p)
                seen.add(key)
            if len(best) >= 2: break
        return best


class MultiplexOptimizer:
    def __init__(self, designer: MGBProbeDesigner):
        self.designer = designer
        self.fluors = [Fluorophore.FAM, Fluorophore.VIC, Fluorophore.NED, Fluorophore.ROX, Fluorophore.CY5]

    def assign_fluorophores_optimal(self, probes: List[MGBProbe]) -> bool:
        used = set()
        i = 0
        while i < len(probes) - 1:
            avail = [f for f in self.fluors if f not in used]
            if len(avail) < 2: return False
            probes[i].fluorophore = avail[0]
            probes[i+1].fluorophore = avail[1]
            used.update(avail[:2])
            i += 2
        return True

    def find_best_4plex(self, all_pairs: List[List[MGBProbe]]) -> Dict:
        best = None
        best_score = -99999
        for combo in combinations(range(len(all_pairs)), min(4, len(all_pairs))):
            selected = []
            for i in combo:
                selected.extend(all_pairs[i][:2])
            if len(selected) >= 4 and self.assign_fluorophores_optimal(selected):
                tms = [p.tm for p in selected]
                score = len(selected) * 100 - (max(tms) - min(tms)) * 10
                if score > best_score:
                    best_score = score
                    best = {
                        "probes": selected,
                        "note": f"Best {len(combo)}-plex found",
                        "tm_range": (min(tms), max(tms))
                    }
        return best or {"error": "No compatible multiplex found"}


def parse_fasta_with_snps(path: str) -> List[Dict]:
    targets = {}
    for rec in SeqIO.parse(path, "fasta"):
        name = rec.id.strip()
        seq = str(rec.seq).upper()
        if "_" in name and name.rsplit("_", 1)[-1] in ("REF", "ALT"):
            base = name.rsplit("_", 1)[0]
            allele = name.split("_")[-1]
            targets.setdefault(base, {})[allele] = seq

    result = []
    for tid, data in targets.items():
        if "REF" in data and "ALT" in data:
            ref = data["REF"]
            alt = data["ALT"]
            for i in range(min(len(ref), len(alt))):
                if ref[i] != alt[i]:
                    result.append({
                        "target_id": tid,
                        "amplicon": ref,
                        "snp_pos": i,
                        "ref": ref[i],
                        "alt": alt[i]
                    })
                    break
    return result


# ------------------- Worker Thread -------------------
class Worker(QThread):
    finished = pyqtSignal(dict)
    error = pyqtSignal(str)
    progress = pyqtSignal(str)

    def __init__(self, path):
        super().__init__()
        self.path = path

    def run(self):
        try:
            self.progress.emit("Parsing FASTA...")
            targets = parse_fasta_with_snps(self.path)
            if not targets:
                self.error.emit("No valid _REF / _ALT pairs found!")
                return

            thermo = ThermodynamicsCalculator()
            designer = MGBProbeDesigner(thermo)
            optimizer = MultiplexOptimizer(designer)
            all_pairs = []

            for t in targets:
                self.progress.emit(f"Designing {t['target_id']}...")
                pair = designer.design_for_snp(t["amplicon"], t["snp_pos"], t["ref"], t["alt"], t["target_id"])
                all_pairs.append(pair)

            result = optimizer.find_best_4plex(all_pairs)
            self.finished.emit(result)

        except Exception as e:
            import traceback
            self.error.emit(f"Error: {str(e)}\n{traceback.format_exc()}")


# ------------------- GUI -------------------
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("MGB Probe Multiplex Designer v2.0")
        self.resize(1300, 800)
        self.result = None
        self.setup_ui()

    def setup_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)

        title = QLabel("MGB Probe Multiplex Designer")
        title.setFont(QFont("Segoe UI", 22, QFont.Weight.Bold))
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(title)

        self.drop = QLabel("Drag & drop FASTA file here\nor click to browse")
        self.drop.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.drop.setStyleSheet("border: 4px dashed #888; padding: 60px; font-size: 18px;")
        self.drop.mousePressEvent = lambda e: self.browse()
        layout.addWidget(self.drop)

        btns = QHBoxLayout()
        self.run_btn = QPushButton("Run Design")
        self.run_btn.clicked.connect(self.start)
        self.run_btn.setEnabled(False)
        self.export_btn = QPushButton("Export to Excel")
        self.export_btn.clicked.connect(self.export)
        self.export_btn.setEnabled(False)
        btns.addWidget(self.run_btn)
        btns.addWidget(self.export_btn)
        layout.addLayout(btns)

        self.status = QLabel("Ready")
        layout.addWidget(self.status)

        self.table = QTableWidget()
        self.table.setColumnCount(9)
        self.table.setHorizontalHeaderLabels([
            "Target", "Allele", "Fluor", "Sequence", "Len", "Tm (°C)", "ΔG", "GC%", "3' Dist"
        ])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        layout.addWidget(self.table)

        self.setAcceptDrops(True)

    def dragEnterEvent(self, e): e.accept() if e.mimeData().hasUrls() else e.ignore()
    def dropEvent(self, e):
        path = e.mimeData().urls()[0].toLocalFile()
        if path.lower().endswith(('.fa', '.fasta', '.txt')):
            self.load_file(path)

    def browse(self):
        path, _ = QFileDialog.getOpenFileName(self, "Open FASTA", "", "FASTA (*.fa *.fasta *.txt)")
        if path: self.load_file(path)

    def load_file(self, path):
        self.fasta_path = path
        self.drop.setText(f"Loaded: {os.path.basename(path)}")
        self.run_btn.setEnabled(True)

    def start(self):
        self.run_btn.setEnabled(False)
        self.status.setText("Processing...")
        self.worker = Worker(self.fasta_path)
        self.worker.finished.connect(self.on_done)
        self.worker.error.connect(lambda m: QMessageBox.critical(self, "Error", m))
        self.worker.progress.connect(self.status.setText)
        self.worker.start()

    def on_done(self, result):
        self.result = result
        self.run_btn.setEnabled(True)
        self.export_btn.setEnabled(True)

        if "error" in result:
            self.status.setText(result["error"])
            return

        probes = result["probes"]
        self.table.setRowCount(len(probes))
        for i, p in enumerate(probes):
            fluor = p.fluorophore.value if p.fluorophore else "-"
            self.table.setItem(i, 0, QTableWidgetItem(p.target_id))
            self.table.setItem(i, 1, QTableWidgetItem(p.target_allele))
            self.table.setItem(i, 2, QTableWidgetItem(fluor))
            self.table.setItem(i, 3, QTableWidgetItem(p.sequence))
            self.table.setItem(i, 4, QTableWidgetItem(str(p.length)))
            self.table.setItem(i, 5, QTableWidgetItem(f"{p.tm:.1f}"))
            self.table.setItem(i, 6, QTableWidgetItem(f"{p.dg:.1f}"))
            self.table.setItem(i, 7, QTableWidgetItem(f"{p.gc_content:.0f}"))
            self.table.setItem(i, 8, QTableWidgetItem(str(p.snp_position)))

        tm_min, tm_max = result["tm_range"]
        self.status.setText(f"{result['note']} | Tm: {tm_min:.1f}–{tm_max:.1f}°C")

    def export(self):
        if not self.result: return
        path, _ = QFileDialog.getSaveFileName(self, "Save", "MGB_Probes.xlsx", "Excel (*.xlsx)")
        if path:
            df = pd.DataFrame([{
                "Target": p.target_id,
                "Allele": p.target_allele,
                "Fluorophore": p.fluorophore.value if p.fluorophore else "",
                "Sequence": p.sequence,
                "Length": p.length,
                "Tm": p.tm,
                "dG": p.dg,
                "GC%": p.gc_content,
                "3'_Dist": p.snp_position
            } for p in self.result["probes"]])
            df.to_excel(path, index=False)
            QMessageBox.information(self, "Saved", f"Exported to {os.path.basename(path)}")


# ------------------- Run -------------------
if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    palette = QPalette()
    palette.setColor(QPalette.ColorRole.Window, QColor(45, 45, 45))
    palette.setColor(QPalette.ColorRole.WindowText, Qt.GlobalColor.white)
    palette.setColor(QPalette.ColorRole.Base, QColor(30, 30, 30))
    palette.setColor(QPalette.ColorRole.Text, Qt.GlobalColor.white)
    palette.setColor(QPalette.ColorRole.Button, QColor(60, 60, 60))
    palette.setColor(QPalette.ColorRole.ButtonText, Qt.GlobalColor.white)
    palette.setColor(QPalette.ColorRole.Highlight, QColor(70, 130, 255))
    app.setPalette(palette)

    win = MainWindow()
    win.show()
    sys.exit(app.exec())
