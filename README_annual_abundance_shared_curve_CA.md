# Abundància amb corba de vol compartida entre anys

`annual_abundance_shared_curve.R` és una anàlisi de sensibilitat separada.
Carregar-lo amb `source()` només defineix funcions: no inicia ajustos ni modifica
els índexs o els models fenològics existents.

## Mètode i límits

Per cada espècie × transsecte, s'ajusta un GAM binomial negatiu amb una forma
estacional comuna i intercepts anuals als anys amb deteccions. Es completen amb
zeros les visites sense registre de l'espècie. Es mantenen els ancoratges hivernals
del procediment original (dies 1–30 i 335–365); no són observacions reals.
La corba es normalitza perquè la suma diària sigui 1. L'índex de cada any és
la suma dels comptatges dividida per la suma de la densitat relativa als dies
visitats. És un estimador corregit per les dates de mostreig, no una suma anual
de comptatges ni una estimació directa de la mida poblacional.

La forma compartida permet incloure anys amb una o dues dates positives. Els
anys sense deteccions obtenen índex zero si la cobertura és adequada; això no
demostra extinció. Cal confirmar que les visites impliquen registre complet
de l'espècie abans d'interpretar els zeros com a no-deteccions vàlides.

Els criteris de cobertura predeterminats, configurables i pendents de sensibilitat,
són: almenys 10 dates de visita; l'interval entre primera i última visita cobreix
almenys el 80% de la massa de la corba; cap buit superior a 28 dies dins del tram
central 5–95%. Els anys que fallen conserven `index_raw`, però `ABUND_INDEX=NA`.
No basta comptar deu visites sense revisar la cobertura estacional.

La corba comuna assumeix forma i calendari constants entre anys. Els desplaçaments
fenològics poden afectar aquesta correcció; cal comparar els índexs nous amb els
originals i inspeccionar les corbes abans d'adoptar-los. No s'estima ni es propaga
la incertesa de la corba o de l'índex. El Gamma demogràfic existent no admet els
zeros nous i requerirà revisió. Aquest script no reajusta models demogràfics.

## Execució explícita

```r
source('annual_abundance_shared_curve.R')
# Pilot petit; cada execució crea una carpeta nova.
pilot <- shared_curve_abundance$run(root='.', max_populations=5)
# Després de revisar el pilot:
# result <- shared_curve_abundance$run(root='.')
```

Per defecte, selecciona poblacions amb almenys deu anys d'índex positiu a
`output/pheno_abund_estimates_allspp.csv`, abans del retall superior. No restringeix
automàticament a les poblacions dels models demogràfics finals. Per fixar aquesta
mostra, passar `population_file` amb columnes `SPECIES,SITE_ID`.
Requereix transsectes amb almenys deu anys amb deu visites i entrada climàtica,
com el flux original. Visites diferents al mateix transsecte-data en anys elegibles
exclouen el transsecte de l'ajust i queden registrades a `failures.csv`: cal resoldre
la unitat de mostreig abans d'ajustar aquestes poblacions.

Els resultats es desen en una carpeta nova sota `output/`: índexs anuals, cobertura,
corbes diàries, correspondència de poblacions, errors, configuració i versions.
No hi ha filtre superior del percentil 99.9. Els anys de traspàs mapen el dia 366
al 365 per compatibilitat amb la graella original.

## Verificació

```sh
Rscript --vanilla tests/test_annual_abundance_shared_curve.R
```

La prova inclou un any amb una detecció, un any sense deteccions, cobertura
insuficient i un buit llarg de visites. No ajusta dades reals.
