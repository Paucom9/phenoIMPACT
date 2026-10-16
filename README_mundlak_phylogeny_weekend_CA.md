# Cinc models Mundlak filogenètics amb malla de 15 km

Scripts preparats el 15 d'octubre de 2026. No contenen dades del projecte i no ajusten cap model en carregar-los.

| Família | Models | Component filogenètic afegit |
| --- | --- | --- |
| Plasticitat–factors ambientals | onset, offset_univoltine, offset_multivoltine | Pendent d'anomalia climàtica per espècie |
| Plasticitat–tendències demogràfiques | univoltine, multivoltine | Pendent temporal per espècie |

Es conserva la formulació Mundlak, els efectes independents d'espècie, la mateixa mostra i la malla de 15 km de cada ajustament independent. Els models demogràfics conserven Gamma/log, AR1 de població amb salts d'any, interceptes i pendents d'espècie i lloc, i camp espacial anual. Els predictors de plasticitat de la primera etapa es mantenen iguals: aquesta comparació avalua l'efecte d'afegir covariància filogenètica al model demogràfic.

## Instal·lació dels scripts

Descomprimeix els cinc fitxers R d'aquest paquet directament a:

`E:/phenoIMPACT project/code/phenoIMPACT`

El paquet inclou:

- `11c_environment_mundlak_phylogeny_15km.R`: tres sensibilitats ambientals.
- `10g_mundlak_demography_phylogeny_15km.R`: dues sensibilitats demogràfiques.
- `run_mundlak_phylogeny_weekend.R`: auditoria, cua i seguiment.
- `pheno_phylogeny_qc_15km_v2.R` i `pheno_plasticity_phylogeny_pilot_15km.R`: dependències de l'11c. La cua les carrega internament de manera que no inicia el pilot antic.

Si R està executant els models independents, deixa que acabin. Els nous scripts escriuen en carpetes de resultats separades. No cal tornar a executar l'11b ni el 10f per carregar aquesta cua.

S'utilitzen els paquets ja instal·lats (`TMB`, `Matrix`, `sdmTMB`, `ape`, `digest`, `callr`, `ps` i dependències). No s'actualitza cap paquet. La primera compilació C++ demogràfica requereix el Rtools corresponent a la versió d'R; per a R 4.5, Rtools45.

## Comprovació inicial: es pot fer abans que acabin els independents

A la consola d'RStudio:

```r
root <- "E:/phenoIMPACT project/code/phenoIMPACT"
source(file.path(root, "run_mundlak_phylogeny_weekend.R"), encoding = "UTF-8")
audit <- audit_mundlak_phylogeny_weekend(project_root = root)
```

Aquesta auditoria només llegeix metadades petites. No obre els models, no compila i no ajusta res. Mentre l'onset independent continuï en marxa és normal que alguns `ready` siguin `FALSE`. La presència d'una carpeta o d'un `prepared_input.rds` no compta com un ajustament complet.

Fonts compatibles:

- Ambientals: `output/phenology_plasticity/environment_mundlak_qc_final_15km/run_*`, versions 1.1/1.2/1.3 i estat final `FIT_OK`. Els offsets han de procedir dels ajustaments post-QC.
- Demogràfiques: `output/population_trends/mundlak_geography_plasticity_qc_final_spatial_15km/run_*`, versió 2.0.2 del 10e/10f i controls numèrics acceptats. Es verifica la mostra FINAL-QC revisada.

Si hi ha més d'una font independent vàlida per a un model, especifica els dos runs explícitament:

```r
audit <- audit_mundlak_phylogeny_weekend(
  project_root = root,
  environment_source_run = "E:/.../environment_mundlak_qc_final_15km/run_...",
  demography_source_run = "E:/.../mundlak_geography_plasticity_qc_final_spatial_15km/run_..."
)
```

## Inici del cap de setmana

Quan els cinc `ready` siguin `TRUE` i els workers independents hagin acabat:

```r
weekend_job <- start_mundlak_phylogeny_weekend(
  project_root = root,
  check_only = FALSE
)
watch_mundlak_phylogeny_weekend(weekend_job)
```

Si has especificat fonts a l'auditoria, passa també `environment_source_run` i `demography_source_run` a l'inici. L'auditoria no modifica els valors per defecte de l'inici.

Ordre: onset ambiental → offset univoltí ambiental → offset multivoltí ambiental → tendències univoltines → tendències multivoltines. Hi ha un controlador lleuger i un sol worker de model alhora. No s'inicia automàticament cap ajustament independent que falti. `start_mundlak_phylogeny_weekend()` sense `check_only = FALSE` només fa l'auditoria.

Mantén RStudio obert i Windows sense suspensió/hibernació. Cada família comprova la RAM disponible abans d'iniciar el seu worker: 30 GiB per als ambientals i 25 GiB per als demogràfics. No és una garantia que la memòria sigui suficient per a tots els ajustaments. El temps i la RAM màxims amb les dades completes encara no s'han mesurat.

## Seguiment, aturada i represa

```r
status_mundlak_phylogeny_weekend(weekend_job)
watch_mundlak_phylogeny_weekend(weekend_job)

# Només si vols aturar expressament la cua i el model actiu:
stop_mundlak_phylogeny_weekend(weekend_job)
```

Interrompre `watch` amb Esc/Ctrl+C pausa el monitor; el worker continua. Tancar RStudio sí que pot aturar els subprocessos.

Després d'una interrupció, torna a carregar el script de la cua i repeteix la mateixa crida `start_mundlak_phylogeny_weekend(..., check_only = FALSE)`, amb els mateixos arguments i fonts. Cada model valida els seus checkpoints: reutilitza els fits complets i recupera la preparació o els paràmetres de l'optimització quan estan disponibles. No es decideix només a partir del resum de la cua.

Es desen la preparació portable, avaluacions finites durant l'optimització, el resultat de l'optimitzador **abans** del càlcul d'incertesa, i el fit final. L'interval de checkpoint de l'optimització és de 120 segons quan una avaluació acaba. Un càlcul natiu llarg, com la construcció de l'Hessià, no pot desar-se a mig executar: després d'una caiguda aquest pas es repeteix.

Per defecte, un error o uns controls numèrics no acceptats queden registrats i la cua continua amb els models següents. Per aturar-se al primer problema, utilitza `stop_on_error = TRUE` en iniciar-la. Els models complets no es tornen a optimitzar només perquè la cua es reprengui; un fit complet però no acceptat necessita revisar-ne el diagnòstic abans de decidir un refinament.

## Resultats i interpretació

- Cua i registres: `output/mundlak_phylogeny_weekend/`.
- Ambientals: `output/phenology_plasticity/environment_mundlak_phylogeny_15km/`.
- Demogràfics: `output/population_trends/mundlak_demography_phylogeny_15km/`.

El monitor mostra el fitxer de registre i la carpeta de resultats del model actiu. `queue_results_latest.csv` recull l'estat dels cinc ajustaments. Cada model exporta controls numèrics, comparacions dels coeficients i errors estàndard, components de variància i comparació ML/AIC. En els demogràfics també es comparen els quatre efectes de plasticitat sobre la tendència (onset/offset × within/between) i els contrastos within–between amb la covariància completa.

Les inferències només s'exporten si passen els controls: convergència, Hessià/covariància definits positius i norma del gradient escalat per covariància inferior a 0.05. També es conserva el gradient brut com a diagnòstic. La variància filogenètica nul·la és una frontera: no s'assigna automàticament un p-valor chi-quadrat estàndard a la seva comparació. L'anàlisi utilitza un únic arbre; no propaga incertesa de l'arbre ni de les plasticitats estimades a la primera etapa.

La validació dels scripts utilitza casos sintètics petits i proves de control de la cua. Els cinc models amb les dades reals del projecte encara s'han d'executar al teu ordinador.

## Validació executada

En R 4.5.3/Linux, TMB 1.9.19 i sdmTMB 1.1.0:

- Versemblança demogràfica: 32 comprovacions numèriques, inclosa equivalència amb el model independent quan es desactiva la filogènia, AR1 amb anys absents i gradients.
- Fonts demogràfiques: 31 comprovacions de signatures, mostra, disseny i rebuig d'entrades inconsistents; exportació estadística: 99 comprovacions.
- Driver demogràfic complet: 24 comprovacions sobre ajustaments sintètics dels dos grups, fit guardat, exports, reutilització en una sessió nova i recuperació del checkpoint anterior al càlcul d'incertesa.
- Cua: 34 comprovacions simulades dels cinc models, errors, caigudes del worker, represa i aturada; 9 comprovacions addicionals amb el controlador i els workers reals de les dues famílies sobre fits sintètics petits.
- Versió ambiental: comprovada amb les sis formulacions onset/offset de les versions independents 1.1/1.2/1.3, ajustaments sintètics i reutilització dels fits.

Aquestes proves no mesuren el cost de les dades completes i no constitueixen una execució a Windows.
