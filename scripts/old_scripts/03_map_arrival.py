#!/usr/bin/env python3
"""步骤3: GMT台站到时分布图"""
import sys, os, yaml, argparse, subprocess, tempfile
for p in ['/home/philiay/shared/codex', '/home/philiay/shared/codex/t_phase_analysis']:
    if p not in sys.path: sys.path.insert(0, p)
from geographiclib.geodesic import Geodesic
import numpy as np

CONFIG = yaml.safe_load(open("/home/sy/codex/config/events.yaml"))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('event')
    args = p.parse_args()

    cfg = CONFIG[args.event]
    case_dir = cfg.get('case_dir', cfg['label'])
    data_dir = f"/mnt/f/sharedate/case/{case_dir}"

    sta_file = f"{data_dir}/gmt_stations.txt"
    if not os.path.exists(sta_file):
        print(f"请先跑 02_extract_stations: {sta_file} 不存在")
        return

    # Read arrivals
    rows = []
    with open(sta_file) as f:
        for line in f:
            if line.startswith('#'): continue
            p = line.split()
            rows.append({'lon': float(p[0]), 'lat': float(p[1]),
                         'tt': float(p[2]), 'sta': p[3]})

    t_min = int(np.floor(min(r['tt'] for r in rows)))
    t_max = int(np.ceil(max(r['tt'] for r in rows)))
    print(f"{cfg['label']}: {len(rows)} stations, {t_min}-{t_max} min")

    # Generate GMT script
    gmt_script = f"""#!/bin/bash
set -e
EVLO={cfg['lon']}; EVLA={cfg['lat']}
REGION="-R-85/-45/20/50"
PROJ="-JM18c"
OUTDIR="{data_dir}"
mkdir -p $OUTDIR; cd $OUTDIR

gmt grdcut @earth_relief_30s $REGION -Gtopo.nc 2>/dev/null || true
gmt makecpt -Cgeo -T-8000/4000/100 -Z > topo.cpt
gmt makecpt -Chot -T{t_min}/{t_max}/0.5 -Z -I > arrival.cpt

gmt set PS_MEDIA 26cx24c PS_PAGE_ORIENTATION portrait

# terrain
gmt grdimage topo.nc $REGION $PROJ -Ctopo.cpt -I+d -t30 -K > map.ps
gmt pscoast $REGION $PROJ -W0.3p,gray40 -Df -A500 -K -O >> map.ps

# stations
awk '!/^#/{{print $1,$2,$3}}' "{sta_file}" | gmt psxy $REGION $PROJ -St0.22c -Carrival.cpt -K -O >> map.ps

# labels
awk '!/^#/{{printf "%.4f %.4f %s\\n", $1+0.3, $2-0.25, $4}}' "{sta_file}" | \\
    gmt pstext $REGION $PROJ -F+f3.5p,Helvetica,gray40+jLM -Gwhite@40 -K -O >> map.ps

# epicenter
echo "$EVLO $EVLA" | gmt psxy $REGION $PROJ -Sa0.5c -Gred -W0.5p,black -K -O >> map.ps
echo "$EVLO $(echo $EVLA+1|bc) Mw {cfg['mw']}" | \\
    gmt pstext $REGION $PROJ -F+f8p,Helvetica-Bold,red -K -O >> map.ps

# frame
gmt psbasemap $REGION $PROJ -Bxa5f1 -Bya5f1 \\
    -BWSne+t"T-wave arrival | {cfg['date']} Mw {cfg['mw']}" -K -O >> map.ps

# colorbar
gmt psscale -Dx5c/0.5c+w10c/0.35c+h+ml -Carrival.cpt -Bxa2+l"Arrival (min)" -O >> map.ps

gmt psconvert map.ps -Tf -A -P; gmt psconvert map.ps -Tg -A -P
echo "Done: map.pdf map.png"
"""
    with open(f"{data_dir}/build_map.sh", 'w') as f:
        f.write(gmt_script)
    subprocess.run(['bash', f"{data_dir}/build_map.sh"])
    print(f"\n✅ {data_dir}/map.png")

if __name__ == '__main__':
    main()
