#!/usr/bin/env python3
"""
MGB Probe Multiplex Designer v2.0
Academic Release — Single file, fully self-contained
GitHub: https://github.com/yashul12-coder/MGB-Probe-Designer
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
            allele = name
