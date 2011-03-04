#!/bin/bash
if [ $# -ne 1 ]; then
	echo "Usage: `basename $0` <file>"
	exit;
fi
perl -p -i -e 's/atcc43049/Haloarcula_marismortui_ATCC_43049/g' $1
perl -p -i -e 's/dsm371/Halobacterium_salinarum_R1/g' $1
perl -p -i -e 's/dsm12286/Halomicrobium_mukohataei_DSM_12286/g' $1
perl -p -i -e 's/dsm16790/Haloquadratum_walsbyi/g' $1
perl -p -i -e 's/dsm12940/Halorhabdus_utahensis_DSM_12940/g' $1
perl -p -i -e 's/atcc49239/Halorubrum_lacusprofundi_ATCC_49239/g' $1
perl -p -i -e 's/atcc33799/Haloarcula_californiae/g' $1
perl -p -i -e 's/atcc33800/Haloarcula_sinaiiensis/g' $1
perl -p -i -e 's/atcc35960/Haloferax_denitrificans/g' $1
perl -p -i -e 's/atcc33500/Haloferax_mediterranei/g' $1
perl -p -i -e 's/atccBAA-1512/Haloferax_mucosum/g' $1
perl -p -i -e 's/atccBAA-897/Haloferax_sulfurifontis/g' $1
perl -p -i -e 's/dsm3757/Haloferax_volcanii_DS2/g' $1
perl -p -i -e 's/jcm10024/Haloarcula_aidinensis_JCM_10024/g' $1
perl -p -i -e 's/jcm12983/Halobiforma_lacisalsi_JCM_12983/g' $1
perl -p -i -e 's/jcm13917/Haloferax_larsenii_JCM_13917/g' $1
perl -p -i -e 's/dsm18796/Halalkalicoccus_jeotgali_B3/g' $1
perl -p -i -e 's/dsm5511/Haloterrigena_turkmenica_DSM_5511/g' $1
perl -p -i -e 's/dsm11551/Halogeometricum_borinquense_DSM_11551/g' $1
perl -p -i -e 's/jcm14265/Halorubrum_ejinorense_JCM_14265/g' $1
perl -p -i -e 's/atcc700873/Halorubrum_hochstenium_ATCC_700873/g' $1
perl -p -i -e 's/atccBAA-652/Haloarcula_sp_GUBF-9_ATCC_BAA_652/g' $1
perl -p -i -e 's/dsm1307/Halococcus_morrhuae_DSM_1307/g' $1
perl -p -i -e 's/dsm8989/Halococcus_salifodinae_DSM_8989/g' $1
perl -p -i -e 's/jcm12892/Halococcus_hamelinensis_JCM_12892/g' $1
perl -p -i -e 's/jcm13561/Halorubrum_litoreum_JCM_13561/g' $1
perl -p -i -e 's/jcm13563/Haloterrigena_limicola_JCM_13563/g' $1
perl -p -i -e 's/jcm13891/Haloterrigena_salina_JCM_13891/g' $1
perl -p -i -e 's/jcm14978/Halorubrum_kocurii_JCM_14978/g' $1
perl -p -i -e 's/jcm10879/Halobiforma_nitratireducens_JCM_10879/g' $1
perl -p -i -e 's/jcm13562/Haloterrigena_longa_JCM_13562/g' $1
perl -p -i -e 's/dsm21966/Halorubrum_xinjiangese_DSM_21966/g' $1
perl -p -i -e 's/jcm13557/Haloarcula_amylolytica/g' $1
perl -p -i -e 's/dsm12282/Haloarcula_argentinensis/g' $1
perl -p -i -e 's/dsm6131/Haloarcula_japonica/g' $1
perl -p -i -e 's/atccBAA-651/Haloarcula_sp_GUBF-8/g' $1
perl -p -i -e 's/jcm8877/Haloarcula_valismortis/g' $1
perl -p -i -e 's/dsm15987/Halobacterium_noricense/g' $1
perl -p -i -e 's/dsm5350/Halococcus_saccharolyticus/g' $1
perl -p -i -e 's/jcm13552/Halococcus_thailandensis/g' $1
perl -p -i -e 's/jcm10717/Haloferax_alexandrinus/g' $1
perl -p -i -e 's/atccBAA-1513/Haloferax_elongans/g' $1
perl -p -i -e 's/atcc33959/Haloferax_gibonsii/g' $1
perl -p -i -e 's/dsm14919/Haloferax_lucentense/g' $1
perl -p -i -e 's/dsm18310/Haloferax_prahovense/g' $1
perl -p -i -e 's/atccBAA-644/Haloferax_sp_GUBF-1/g' $1
perl -p -i -e 's/atccBAA-645/Haloferax_sp_GUBF-2/g' $1
perl -p -i -e 's/atccBAA-646/Haloferax_sp_GUBF-3/g' $1
perl -p -i -e 's/jcm13560/Halorubrum_aidingense/g' $1
perl -p -i -e 's/jcm13916/Halorubrum_arcis/g' $1
perl -p -i -e 's/dsm19288/Halorubrum_californiensis/g' $1
perl -p -i -e 's/dsm10284/Halorubrum_coriense/g' $1
perl -p -i -e 's/jcm10118/Halorubrum_distributum/g' $1
perl -p -i -e 's/jcm9100/Halorubrum_distributum/g' $1
perl -p -i -e 's/dsm21995/Halorubrum_lipolyticum/g' $1
perl -p -i -e 's/dsm1137/Halorubrum_saccharovorum/g' $1
perl -p -i -e 's/dsm14210/Halorubrum_tebenquichense/g' $1
perl -p -i -e 's/jcm10247/Halorubrum_terrestre/g' $1
perl -p -i -e 's/jcm14848/Halosarcina_pallida/g' $1
perl -p -i -e 's/jcm11222/Halosimplex_carlsbadense/g' $1
perl -p -i -e 's/dsm11522/Haloterrigena_thermotolerans/g' $1
perl -p -i -e 's/jcm14624/Halovivax_asiaticus/g' $1
perl -p -i -e 's/dsm2160/Natronomonas_pharaonis_DSM_2160/g' $1
perl -p -i -e 's/dsm12281/Natrialba_taiwanensis_DSM_12281/g' $1
perl -p -i -e 's/atcc43099/Natrialba_magadii_ATCC_43099/g' $1
perl -p -i -e 's/dsm10524/Natronococcus_amylolyticus_DSM_10524/g' $1
perl -p -i -e 's/dsm13077/Natrialba_aegyptia_DSM_13077/g' $1
perl -p -i -e 's/jcm14089/Natronorubrum_sulfidifaciens_JCM_14089/g' $1
perl -p -i -e 's/dsm3751/Natrinema_pallidum_DSM_3751/g' $1
perl -p -i -e 's/dsm12278/Natrialba_asiatica/g' $1
perl -p -i -e 's/jcm10990/Natrialba_chahannoensis/g' $1
perl -p -i -e 's/jcm10989/Natrialba_hulunbeirensis/g' $1
perl -p -i -e 's/dsm3394/Natrialba_magadii/g' $1
perl -p -i -e 's/jcm12890/Natrinema_altunense/g' $1
perl -p -i -e 's/jcm13890/Natrinema_ejinorense/g' $1
perl -p -i -e 's/jcm14663/Natrinema_gari/g' $1
perl -p -i -e 's/dsm15624/Natrinema_pellirubrum/g' $1
perl -p -i -e 's/jcm10478/Natrinema_versiforme/g' $1
perl -p -i -e 's/dsm3393/Natronobacterium_gregoryi/g' $1
perl -p -i -e 's/dsm18795/Natronococcus_jeotgali/g' $1
perl -p -i -e 's/jcm12255/Natronolimnobius_innermongolicus/g' $1
perl -p -i -e 's/jcm10635/Natronorubrum_bangense/g' $1
perl -p -i -e 's/jcm10636/Natronorubrum_tibetense/g' $1
perl -p -i -e 's/jcm11081/Halobacterium_sp_NRC-1/g' $1
