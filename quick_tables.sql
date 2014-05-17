-- psql -d IrisQuebec -f quick_tables.sql
-- SQL
SET client_min_messages TO 'warning';

drop table if exists aman.chandan_drugs;
create table aman.chandan_drugs as 
select x.id,x.discharge, 'drug' || y.cod_denom_comne || '_' || z.nom 
from aman.readmissions_top20 x
left join admission_drugs y
on x.id=y.id and x.discharge=y.discharge
left join metadata.denominations_communes z
on y.cod_denom_comne=z.code
group by x.id,x.discharge,y.cod_denom_comne,z.nom
order by x.id,x.discharge,y.cod_denom_comne,z.nom;

drop table if exists aman.chandan_diagnoses;
create table aman.chandan_diagnoses as 
select x.id,x.discharge,'diag' || y.icd9 || '_' || COALESCE(a.icd9_desc,i.nom) as icd9
from aman.readmissions_top20 x
left join tight_icd9 y
on x.id=y.id and x.discharge=y.discharge
left join metadata.classification_internationale_des_maladies_9e_revision_quebec i
on y.icd9=i.code
left join aman.icd9_american a on i.code=a.icd9
group by x.id,x.discharge,y.icd9,a.icd9_desc,i.nom
order by x.id,x.discharge,y.icd9,a.icd9_desc,i.nom;

drop table if exists aman.chandan_procedures;
create table aman.chandan_procedures as 
select x.id,x.discharge,'proc' || y.ccp || '_' || y.ccp_desc
from aman.readmissions_top20 x
left join id_discharge_ccp_full y
on x.id=y.id and x.discharge=y.discharge
group by x.id,x.discharge,y.ccp,y.ccp_desc
order by x.id,x.discharge,y.ccp,y.ccp_desc;


/*
drop table if exists aman.chandan_drugs_ahfs;
create table aman.chandan_drugs_ahfs as 
select x.id,x.discharge,z.ahfs_class_descr
                   from aman.readmissions_top20 x
                   left join admission_drugs y
                   on x.id=y.id and x.discharge=y.discharge
                   left join drugs.din2class z
                   on y.cod_din=z.din
                   group by x.id,x.discharge,z.ahfs_class_descr;
*/