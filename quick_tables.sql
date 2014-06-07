-- psql -d IrisQuebec -f quick_tables.sql
-- SQL
SET client_min_messages TO 'warning';

drop table if exists aman.drugs_one_year;
create table aman.drugs_one_year as
select x.id, x.discharge, 'drug' || y.cod_denom_comne || '_' || z.nom  as drug
from aman.readmissions_top20_one_year x join data.ramq_i_iqa2_dem_paimt_med y on 
x.id=y.no_iciq_ben AND 
((y.dat_serv between (x.admit-365) AND x.admit) OR
((y.dat_serv+y.nbr_jr_duree_trait) between (x.admit-365) AND x.admit))
left join metadata.denominations_communes z
on y.cod_denom_comne=z.code
group by x.id,x.discharge,y.cod_denom_comne,z.nom;

drop table if exists aman.diags_one_year;
create table aman.diags_one_year as 
select x.id,x.discharge,'diag' || y.icd9 || '_' || COALESCE(a.icd9_desc,i.nom) as icd9
from aman.readmissions_top20_one_year x
left join tight_icd9 y
on x.id=y.id and (x.discharge>y.discharge AND (y.discharge>=x.discharge-365))
left join metadata.classification_internationale_des_maladies_9e_revision_quebec i
on y.icd9=i.code
left join aman.icd9_american a on i.code=a.icd9
group by x.id,x.discharge,y.icd9,a.icd9_desc,i.nom
order by x.id,x.discharge,y.icd9,a.icd9_desc,i.nom;

drop table if exists aman.procs_one_year;
create table aman.procs_one_year as 
select x.id,x.discharge,'proc' || y.ccp || '_' || y.ccp_desc
from aman.readmissions_top20_one_year x
left join id_discharge_ccp_full y
on x.id=y.id and (x.discharge>y.discharge AND (y.discharge>=x.discharge-365))
group by x.id,x.discharge,y.ccp,y.ccp_desc
order by x.id,x.discharge,y.ccp,y.ccp_desc;