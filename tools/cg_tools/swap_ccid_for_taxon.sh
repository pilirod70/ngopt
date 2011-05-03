#!/bin/bash
if [ $# -ne 1 ]; then
    echo "Usage: `basename $0` <file>"
    exit;
fi

perl -p -i -e "s/dsm18195/367189/g" $1 
perl -p -i -e "s/dsm18796/413810/g" $1 
perl -p -i -e "s/jcm11890/175632/g" $1 
perl -p -i -e "s/jcm10024/56545/g" $1 
perl -p -i -e "s/jcm13557/396317/g" $1 
perl -p -i -e "s/dsm12282/43776/g" $1 
perl -p -i -e "s/jcm8911/jcm8911/g" $1 
perl -p -i -e "s/dsm6131/29282/g" $1 
perl -p -i -e "s/dsm12286/dsm12286/g" $1 
perl -p -i -e "s/atcc33799/662475/g" $1 
perl -p -i -e "s/atcc33800/662476/g" $1 
perl -p -i -e "s/atccBAA[-,_]651/atccBAA-651/g" $1 
perl -p -i -e "s/atccBAA[-,_]652/atccBAA-652/g" $1 
perl -p -i -e "s/jcm8877/28442/g" $1 
perl -p -i -e "s/dsm15987/223182/g" $1 
perl -p -i -e "s/dsm9297/43928/g" $1 
perl -p -i -e "s/jcm11627/148448/g" $1 
perl -p -i -e "s/jcm12983/229731/g" $1 
perl -p -i -e "s/jcm10879/130048/g" $1 
perl -p -i -e "s/dsm14522/179637/g" $1 
perl -p -i -e "s/jcm12892/332168/g" $1 
perl -p -i -e "s/dsm13077/129789/g" $1 
perl -p -i -e "s/dsm1307/931277/g" $1 
perl -p -i -e "s/jcm13587/224402/g" $1 
perl -p -i -e "s/dsm5350/62319/g" $1 
perl -p -i -e "s/dsm8989/36738/g" $1 
perl -p -i -e "s/jcm13552/335952/g" $1 
perl -p -i -e "s/jcm10717/114529/g" $1 
perl -p -i -e "s/atcc35960/662478/g" $1 
perl -p -i -e "s/atccBAA[-,_]1513/403191/g" $1 
perl -p -i -e "s/atcc33959/35746/g" $1 
perl -p -i -e "s/jcm13917/302484/g" $1 
perl -p -i -e "s/dsm14919/2254/g" $1 
perl -p -i -e "s/atcc33500/523841/g" $1 
perl -p -i -e "s/atccBAA[-,_]1512/662479/g" $1 
perl -p -i -e "s/dsm18310/381852/g" $1 
perl -p -i -e "s/atcc51408/atcc51408/g" $1 
perl -p -i -e "s/atccBAA[-,_]644/atccBAA-644/g" $1 
perl -p -i -e "s/atccBAA[-,_]645/atccBAA-645/g" $1 
perl -p -i -e "s/atccBAA[-,_]646/atccBAA-646/g" $1 
perl -p -i -e "s/atccBAA[-,_]897/662480/g" $1 
perl -p -i -e "s/dsm3757/309800/g" $1 
perl -p -i -e "s/dsm11551/469382/g" $1 
perl -p -i -e "s/jcm14033/387343/g" $1 
perl -p -i -e "s/dsm17983/376171/g" $1 
perl -p -i -e "s/jcm13560/368623/g" $1 
perl -p -i -e "s/jcm12358/261290/g" $1 
perl -p -i -e "s/jcm13916/368454/g" $1 
perl -p -i -e "s/dsm19288/416585/g" $1 
perl -p -i -e "s/dsm10284/64713/g" $1 
perl -p -i -e "s/jcm9100/29283/g" $1 
perl -p -i -e "s/jcm10118/29283/g" $1 
perl -p -i -e "s/jcm14265/425309/g" $1 
perl -p -i -e "s/atcc700873/atcc700873/g" $1 
perl -p -i -e "s/jcm14978/478441/g" $1 
perl -p -i -e "s/dsm21995/368624/g" $1 
perl -p -i -e "s/jcm13561/409316/g" $1 
perl -p -i -e "s/dsm1137/2248/g" $1 
perl -p -i -e "s/dsm3755/dsm3755/g" $1 
perl -p -i -e "s/dsm14210/119434/g" $1 
perl -p -i -e "s/jcm10247/267368/g" $1 
perl -p -i -e "s/jcm11888/175631/g" $1 
perl -p -i -e "s/dsm8800/63740/g" $1 
perl -p -i -e "s/dsm21966/261291/g" $1 
perl -p -i -e "s/jcm14848/411361/g" $1 
perl -p -i -e "s/jcm11222/171164/g" $1 
perl -p -i -e "s/jcm13563/370323/g" $1 
perl -p -i -e "s/jcm13562/370324/g" $1 
perl -p -i -e "s/jcm12889/301967/g" $1 
perl -p -i -e "s/jcm13891/504937/g" $1 
perl -p -i -e "s/dsm11522/121872/g" $1 
perl -p -i -e "s/dsm5511/543526/g" $1 
perl -p -i -e "s/jcm14624/332953/g" $1 
perl -p -i -e "s/dsm18193/387341/g" $1 
perl -p -i -e "s/dsm_12278/64602/g" $1 
perl -p -i -e "s/jcm10990/68911/g" $1 
perl -p -i -e "s/jcm10989/123783/g" $1 
perl -p -i -e "s/dsm3394/13769/g" $1 
perl -p -i -e "s/dsm12281/160846/g" $1 
perl -p -i -e "s/jcm12890/222984/g" $1 
perl -p -i -e "s/jcm13890/373386/g" $1 
perl -p -i -e "s/jcm14663/419186/g" $1 
perl -p -i -e "s/dsm3751/69527/g" $1 
perl -p -i -e "s/dsm15624/797303/g" $1 
perl -p -i -e "s/jcm10478/88724/g" $1 
perl -p -i -e "s/dsm3393/44930/g" $1 
perl -p -i -e "s/dsm10524/44470/g" $1 
perl -p -i -e "s/dsm18795/413812/g" $1 
perl -p -i -e "s/dsm3396/694430/g" $1 
perl -p -i -e "s/jcm12253/253108/g" $1 
perl -p -i -e "s/jcm12255/253107/g" $1 
perl -p -i -e "s/jcm13488/348826/g" $1 
perl -p -i -e "s/jcm10635/61858/g" $1 
perl -p -i -e "s/jcm14089/388259/g" $1 
perl -p -i -e "s/jcm10636/63128/g" $1 
