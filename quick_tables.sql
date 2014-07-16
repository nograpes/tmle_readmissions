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

-- select * from aman.icd9_outpatient_codes limit 100; -- id, date, icd9
-- 6.1 million
drop table if exists temp_outpatient_diags_one_year;
create table temp_outpatient_diags_one_year as
select x.id,x.discharge,'diag' || y.icd9 || '_' || COALESCE(a.icd9_american_descr,i.nom) as icd9,a.icd9_american
from aman.readmissions_top20_one_year x
left join icd9_outpatient_codes y
on x.id=y.id 
left join metadata.classification_internationale_des_maladies_9e_revision_quebec i
on y.icd9=i.code
left join aman.icd9_to_icd9_american a on i.code=a.icd9
-- left join aman.icd9_american a on i.code=a.icd9
where (x.discharge>y.date AND (y.date>=x.discharge-365))
group by x.id,x.discharge,y.icd9,a.icd9_american_descr,i.nom,a.icd9_american
order by x.id,x.discharge,y.icd9,a.icd9_american_descr,i.nom,a.icd9_american;


drop table if exists temp_diags_one_year;
create table temp_diags_one_year as 
select x.id,x.discharge,'diag' || y.icd9 || '_' || COALESCE(a.icd9_american_descr,i.nom) as icd9,a.icd9_american
from aman.readmissions_top20_one_year x
left join tight_icd9 y
on x.id=y.id and (x.discharge>y.discharge AND (y.discharge>=x.discharge-365))
left join metadata.classification_internationale_des_maladies_9e_revision_quebec i
on y.icd9=i.code
left join aman.icd9_to_icd9_american a on i.code=a.icd9
-- left join aman.icd9_american a on i.code=a.icd9
group by x.id,x.discharge,y.icd9,a.icd9_american_descr,i.nom,a.icd9_american
order by x.id,x.discharge,y.icd9,a.icd9_american_descr,i.nom,a.icd9_american;

-- select count(*) from temp_diags_one_year; -- 2,854,965
-- select count(*) from temp_outpatient_diags_one_year; -- 6,180,368

drop table if exists aman.diags_one_year;
create table aman.diags_one_year as
select id,discharge,icd9 
from temp_diags_one_year;

drop table if exists aman.diags_w_outpatient_one_year;
create table aman.diags_w_outpatient_one_year as
select * from 
(select * from aman.diags_one_year
UNION ALL
select id,discharge,icd9 from temp_outpatient_diags_one_year)_
group by id,discharge,icd9
order by id,discharge,icd9
;

drop table if exists aman.charlson_comorbidities;
create table aman.charlson_comorbidities (cc text, icd9_american text);
copy aman.charlson_comorbidities from '/home/aman/repo/thesis/data/charlson/charlson.csv' WITH CSV HEADER;

drop table if exists temp_charlson_comorbidities_one_year;
create temp table temp_charlson_comorbidities_one_year as
select x.id,x.discharge,'cc_'||y.cc as cc from
(select id,discharge,icd9_american from temp_diags_one_year
UNION ALL
select id,discharge,icd9_american from temp_outpatient_diags_one_year) x
join aman.charlson_comorbidities y
on x.icd9_american = y.icd9_american
group by id,discharge,cc
order by id,discharge,cc;

drop table if exists aman.charlson_comorbidities_one_year;
create table aman.charlson_comorbidities_one_year as
select x.id,x.discharge,y.cc from aman.readmissions_top20_one_year x
left join temp_charlson_comorbidities_one_year y 
on x.id=y.id and x.discharge=y.discharge;

drop table if exists aman.charlson_comorbidities_counts_one_year;
create table aman.charlson_comorbidities_counts_one_year as
select x.id,x.discharge,y.count from readmissions_top20_one_year x
left join (select id, discharge, count(*) as count from temp_charlson_comorbidities_one_year group by id,discharge) y 
on x.id=y.id and x.discharge=y.discharge;


/*
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



drop table if exists aman.outpatient_diags_one_year;
create table aman.outpatient_diags_one_year as
select x.id,x.discharge,'diag' || y.icd9 || '_' || COALESCE(a.icd9_desc,i.nom) as icd9
from aman.readmissions_top20_one_year x
left join icd9_outpatient_codes y
on x.id=y.id 
left join metadata.classification_internationale_des_maladies_9e_revision_quebec i
on y.icd9=i.code
left join aman.icd9_american a on i.code=a.icd9
where (x.discharge>y.date AND (y.date>=x.discharge-365))
group by x.id,x.discharge,y.icd9,a.icd9_desc,i.nom
order by x.id,x.discharge,y.icd9,a.icd9_desc,i.nom;


select * from aman.outpatient_diags_one_year limit 100;
select * from aman.diags_one_year limit 100;
*/

drop table if exists aman.procs_one_year;
create table aman.procs_one_year as 
select x.id,x.discharge,'proc' || y.ccp || '_' || y.ccp_desc
from aman.readmissions_top20_one_year x
left join id_discharge_ccp_full y
on x.id=y.id and (x.discharge>y.discharge AND (y.discharge>=x.discharge-365))
group by x.id,x.discharge,y.ccp,y.ccp_desc
order by x.id,x.discharge,y.ccp,y.ccp_desc;