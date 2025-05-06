#!/usr/bin/env python
import json
for srm in [0.01,0.1,1,10,100]:
    for nvm in [0.01,0.1,1,10,100]:
        for energy in ["pair zscore","sRNA zscore"]:
            for marker in ["16S-rRNA","rpoB"]:
                config = json.load(open("template.0425.json"))
                config['marker'] = marker
                config['sRNA-ids'] = "sRNA-ids.txt"
                config["tag"] = "202504"
                energy_str = energy.replace(" ",".")
                config['denoise']['srm'] = srm
                config['denoise']['nvm'] = nvm
                if energy == "pair zscore":            
                    config['normalize'] = 1
                config['weights'] = {energy: 0.7, 'hfq Z': 0.3}
                fout = open(f"grids.20250425/{srm}-{nvm}-{energy_str}-{marker}-0.7.json","w")
                json.dump(config,fout)
                fout.close()
                config['weights'] = {energy: 1.0, 'hfq Z': 0.0}
                fout = open(f"grids.20250425/{srm}-{nvm}-{energy_str}-{marker}-1.0.json","w")
                json.dump(config,fout)
                fout.close()
