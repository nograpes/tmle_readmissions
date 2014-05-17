-- psql -d IrisQuebec -f quick_tables.sql
-- SQL
SET client_min_messages TO 'warning';
drop table if exists aman.chandan_drugs;
create table aman.chandan_drugs as 
select x.id,x.discharge,z.nom
                  from aman.readmissions_top20 x
                  left join admission_drugs y
                  on x.id=y.id and x.discharge=y.discharge
                  left join metadata.denominations_communes z
                  on y.cod_denom_comne=z.code;


drop table if exists aman.chandan_diagnoses;
create table aman.chandan_diagnoses as 
select x.id,x.discharge,y.icd9
                      from aman.readmissions_top20 x
                      left join tight_icd9 y
                      on x.id=y.id and x.discharge=y.discharge
                      order by id,discharge;


drop table if exists aman.chandan_procedures;
create table aman.chandan_procedures as 
select x.id,x.discharge,y.ccp_desc
                       from aman.readmissions_top20 x
                       left join id_discharge_ccp_full y
                       on x.id=y.id and x.discharge=y.discharge
                       order by id,discharge;

drop table if exists aman.chandan_drugs_ahfs;
create table aman.chandan_drugs_ahfs as 
select x.id,x.discharge,z.ahfs_class_descr
                   from aman.readmissions_top20 x
                   left join admission_drugs y
                   on x.id=y.id and x.discharge=y.discharge
                   left join drugs.din2class z
                   on y.cod_din=z.din
                   group by x.id,x.discharge,z.ahfs_class_descr;