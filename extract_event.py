"""Extract a single event from a MicroBooNE OpenSamples HDF5 file to JSON.

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
        hits = f["hit_table"]
        a, b = row_range(hits, event_index)
        plane = hits["local_plane"][a:b, 0]
        wire = hits["local_wire"][a:b, 0]
        time = hits["local_time"][a:b, 0]
        charge = hits["integral"][a:b, 0]
        planes = {}
        for p in (0, 1, 2):
            m = plane == p
            planes[str(p)] = {
                "wire": wire[m].tolist(),
                "time": np.round(time[m], 1).tolist(),
                "charge": np.round(charge[m], 1).tolist(),
            }
        out = {"meta": meta, "planes": planes,
               "acknowledgment": "MicroBooNE Collaboration public data sets, "
                                 "DOI 10.5281/zenodo.7261921 (CC-BY 4.0)"}
    with open(out_path, "w") as g:
        json.dump(out, g)
    n = sum(len(planes[p]["wire"]) for p in planes)
    print(f"Wrote {out_path}: {n} hits, E_nu={meta['nu_energy_gev']:.2f} GeV")


if __name__ == "__main__":
    h5 = sys.argv[1]
    idx = int(sys.argv[2])
    out = sys.argv[3] if len(sys.argv) > 3 else f"event_{idx}.json"
    extract(h5, idx, out)
