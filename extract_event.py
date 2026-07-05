"""Extract a single event from a MicroBooNE OpenSamples HDF5 file to JSON.

Includes per-hit truth labels (edep_table join) and true particle
trajectories (particle_table) for 2D truth coloring and 3D display.

Usage: python extract_event.py <file.h5> <event_index> [out.json]
Data: MicroBooNE Public Datasets (CC-BY 4.0), DOI 10.5281/zenodo.7261921
"""
import sys, json
import numpy as np
import h5py


def row_range(group, event_index):
    """Rows for event_index using event_id.seq_cnt [seq, count] contiguous layout."""
    seq = group["event_id.seq_cnt"][:]
    counts = seq[:, 1]
    starts = np.concatenate([[0], np.cumsum(counts)[:-1]])
    sel = np.where(seq[:, 0] == event_index)[0]
    if len(sel) == 0:
        return 0, 0
    i = sel[0]
    return int(starts[i]), int(starts[i] + counts[i])


def extract(path, event_index, out_path):
    with h5py.File(path, "r") as f:
        ev = f["event_table"]
        meta = {
            "source_file": path.split("\\")[-1],
            "event_index": event_index,
            "event_id": [int(x) for x in ev["event_id"][event_index]],
            "nu_pdg": int(ev["nu_pdg"][event_index, 0]),
            "nu_energy_gev": float(ev["nu_energy"][event_index, 0]),
            "lep_energy_gev": float(ev["lep_energy"][event_index, 0]),
            "is_cc": bool(ev["is_cc"][event_index, 0]),
            "nu_vtx_cm": [float(x) for x in ev["nu_vtx"][event_index]],
        }
        # --- truth: map hit_id -> g4_id via edep_table (best energy match) ---
        ed = f["edep_table"]
        ea, eb = row_range(ed, event_index)
        ed_hit = ed["hit_id"][ea:eb, 0]
        ed_g4 = ed["g4_id"][ea:eb, 0]
        ed_en = ed["energy"][ea:eb, 0]
        order = np.argsort(ed_en)  # ascending: later writes win = max energy
        hit2g4 = {}
        for j in order:
            hit2g4[int(ed_hit[j])] = int(ed_g4[j])
        # --- true particles for this event ---
        pt = f["particle_table"]
        pa, pb = row_range(pt, event_index)
        p_g4 = pt["g4_id"][pa:pb, 0]
        p_pdg = pt["g4_pdg"][pa:pb, 0]
        p_start = pt["start_position_corr"][pa:pb]
        p_end = pt["end_position_corr"][pa:pb]
        p_mom = pt["momentum"][pa:pb, 0]
        # keep particles that actually produced hits, plus any with momentum > 50 MeV
        hitmakers = set(hit2g4.values())
        g4_to_slot = {}
        particles = []
        for j in range(len(p_g4)):
            gid = int(p_g4[j])
            if gid not in hitmakers and p_mom[j] < 0.05:
                continue
            g4_to_slot[gid] = len(particles)
            particles.append({
                "g4_id": gid, "pdg": int(p_pdg[j]),
                "mom_gev": round(float(p_mom[j]), 3),
                "start": [round(float(x), 1) for x in p_start[j]],
                "end": [round(float(x), 1) for x in p_end[j]],
                "n_hits": 0,
            })
        # --- hits per plane with truth slot (-1 = unmatched / cosmic data) ---
        hits = f["hit_table"]
        a, b = row_range(hits, event_index)
        plane = hits["local_plane"][a:b, 0]
        wire = hits["local_wire"][a:b, 0]
        time = hits["local_time"][a:b, 0]
        charge = hits["integral"][a:b, 0]
        hid = hits["hit_id"][a:b, 0]
        planes = {}
        for p in (0, 1, 2):
            m = plane == p
            slots = []
            for h in hid[m]:
                s = g4_to_slot.get(hit2g4.get(int(h), -1), -1)
                if s >= 0:
                    particles[s]["n_hits"] += 1
                slots.append(s)
            planes[str(p)] = {
                "wire": wire[m].tolist(),
                "time": np.round(time[m], 1).tolist(),
                "charge": np.round(charge[m], 1).tolist(),
                "truth": slots,
            }
        out = {"meta": meta, "planes": planes, "particles": particles,
               "acknowledgment": "MicroBooNE Collaboration public data sets, "
                                 "DOI 10.5281/zenodo.7261921 (CC-BY 4.0)"}
    with open(out_path, "w") as g:
        json.dump(out, g)
    n = sum(len(planes[p]["wire"]) for p in planes)
    nt = sum(1 for p in planes for s in planes[p]["truth"] if s >= 0)
    print(f"Wrote {out_path}: {n} hits ({nt} truth-matched), "
          f"{len(particles)} particles, E_nu={meta['nu_energy_gev']:.2f} GeV")


if __name__ == "__main__":
    h5 = sys.argv[1]
    idx = int(sys.argv[2])
    out = sys.argv[3] if len(sys.argv) > 3 else f"event_{idx}.json"
    extract(h5, idx, out)
