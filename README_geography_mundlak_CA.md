# Patrons geogràfics de plasticitat: quatre models directes

El fitxer `geography_mundlak_by_voltinism_qc_final.R` prepara quatre ajustos:

| Ajust | Resposta anual | Finestra d’anomalia | Covariable additiva |
|---|---|---|---|
| offset multivoltines | OFFSET_mean | 90 dies | ONSET_mean_z |
| onset multivoltines | ONSET_mean | 60 dies | — |
| offset univoltines | OFFSET_mean | 90 dies | ONSET_mean_z |
| onset univoltines | ONSET_mean | 60 dies | — |

La fórmula fixa de cada ajust és:

```r
fenologia ~ anomalia * (
  latitude_within_species_10deg + latitude_species_mean_10deg +
  longitude_within_species_10deg + longitude_species_mean_10deg
)
# + ONSET_mean_z als offsets
```

Es conserven els termes aleatoris de la fórmula d’origen, la família gaussiana
amb enllaç identitat, ML, el camp espacial anual IID Matérn i la mateixa malla
de 15 km. `spatial = "off"` exclou el camp espacial permanent; el camp anual
es manté amb `spatiotemporal = "iid"` i `time = "YEAR"`.

La resposta és la fenologia observada anual. Els coeficients ambientals i les
plasticitats predites pels models ambientals no intervenen en aquests ajustos.
La configuració ambiental acordada continua sense latitud/longitud explícites.

## Ús

Desa el fitxer R a `E:/phenoIMPACT project/code/phenoIMPACT`.
Per carregar les funcions i comprovar l’estat dels inputs:

```r
source("geography_mundlak_by_voltinism_qc_final.R", encoding = "UTF-8")
pheno_geo_mundlak$plan()
pheno_geo_mundlak$preflight()
```

`source()` només defineix `pheno_geo_mundlak`. `preflight()` llegeix metadades
petites, informa de RAM/disc i processos coneguts, i no carrega cap ajust ni
entrada de dades gran. Aquesta comprovació es pot fer mentre corren els ambientals.

Quan els quatre ambientals siguin `FIT_OK` i el seu procés hagi acabat:

```r
geo_job <- pheno_geo_mundlak$run()
pheno_geo_mundlak$status()
```

El llançador exigeix 35 GiB de RAM disponibles i, per defecte, 50 GiB de disc.
La disponibilitat es torna a comprovar abans de cada ajust. Aquests llindars són
una precaució operativa; no són una mesura del consum màxim dels models.
Es rebutja el llançament si es detecta un altre workflow phenoIMPACT actiu.

El controlador continua en segon pla si es tanca la sessió R que l’ha iniciat.
Cada ajust s’executa en un procés nou, amb un fil, per alliberar la memòria abans
del següent. Després de reobrir R, torna a fer `source(...)` i
`pheno_geo_mundlak$status()` per recuperar l’estat des del disc.

## Inputs i comparabilitat

La font predeterminada és:

```text
output/phenology_plasticity/environment_mundlak_by_voltinism_qc_final_15km/run_de0890c85d4a
```

Només s’accepta la versió ambiental
`environment_mundlak_by_voltinism_qc_final_15km_2.1_resilient`.
Es llegeixen els `prepared_input.rds` dels quatre ajustos completats; no els
seus `fit_mundlak.rds`. Això conserva exactament la mostra després del FINAL-QC
dels offsets i de la classificació de voltinisme dels onsets, sense tornar a
carregar els ajustos originals grans.

El codi comprova la correspondència entre resposta, finestra, grup, covariable
d’onset, ordre de files i projecció de la malla. No elimina files per recuperar
coordenades ni reselecciona finestres. La signatura de cada execució inclou les
empremtes dels inputs, el codi i les versions dels paquets.

Si el mateix workflow ambiental s’executa en una altra carpeta, es pot indicar:

```r
pheno_geo_mundlak$preflight(environment_run = "E:/.../run_correcte")
geo_job <- pheno_geo_mundlak$run(environment_run = "E:/.../run_correcte")
```

## Escales geogràfiques i interpretació

Les coordenades originals en metres EPSG:3035 es comproven contra `x_km/y_km`
i es transformen a graus EPSG:4326. Per a cada eix:

- El component **within** és la coordenada de la població menys la mitjana de
  la seva espècie, dividida per 10.
- El component **between** és la mitjana de l’espècie menys una referència
  comuna, dividida per 10. La referència és la mitjana no ponderada de les
  mitjanes de les espècies; només fixa l’origen de l’escala.
- Cada espècie es resumeix donant el mateix pes a cada població espècie × lloc
  del grup de voltinisme corresponent, independentment dels anys de seguiment.

Les mitjanes representen la distribució **mostrejada** de cada espècie dins el
grup, no el centre de tota la seva distribució geogràfica. Les espècies amb una
sola població no aporten contrast geogràfic dins espècie.

Les interaccions són canvis en la sensibilitat fenològica per 10 graus, ajustats
per l’altre eix i l’altra escala. La unitat de l’anomalia tèrmica es conserva
exactament com a la font; no es reetiqueta automàticament com °C.

L’orientació crua és dies: negatiu indica més avançament i positiu més retard.
També s’exporta una orientació biològica: positiu indica més avançament als
onsets i a l’offset univoltí, i més retard de l’offset multivoltí condicionat a
l’onset anual.

## Resultats i recuperació

Els resultats es desen en una carpeta nova i separada:

```text
output/phenology_plasticity/geography_mundlak_by_voltinism_qc_final_15km/run_<signatura>
```

Fitxers principals:

| Fitxer | Contingut |
|---|---|
| `ALL_geography_within_between_effects.csv` | Quatre gradients per ajust, amb IC i tests Wald |
| `ALL_within_minus_between_contrasts.csv` | Contrast directe entre les dues escales de cada eix |
| `ALL_joint_geography_tests.csv` | Tests Wald conjunts within (2 gl), between (2 gl), total (4 gl) |
| `ALL_fit_checks.csv` | Convergència, Hessiana i gradient |
| `ALL_preparation_checks.csv` | Mostra, rang de la matriu i equivalència dels inputs |
| `geography_mundlak_results_to_review.zip` | Taules i diagnòstics textuals, si està instal·lat `zip` |

Cada subcarpeta inclou la covariància fixa completa, suport geogràfic per espècie,
correlacions i VIF, entrades preparades i `fit_geography.rds`.
`population_fixed_geographic_sensitivity.csv` conté la sensibilitat del component
geogràfic fix als contextos observats, amb incertesa propagada amb tota la
covariància fixa. Exclou els pendents aleatoris d’espècie. Aquestes prediccions
no són observacions independents per fer una altra regressió geogràfica.

En tornar a executar `run()` amb la mateixa configuració es reutilitzen els
ajustos completats. Si falta alguna exportació, es recupera des de l’ajust
guardat. Si el procés s’interromp **durant l’optimització abans de guardar l’ajust**,
aquell ajust torna a començar amb els valors inicials; no hi ha recuperació
de la darrera iteració. Els ajustos anteriors completats es conserven.

`FIT_OK` exigeix codi de convergència 0, Hessiana definida positiva i gradient
màxim absolut inferior a 0,001. Només indica comprovacions numèriques satisfactòries.
La revisió dels residus i de la dependència espacial residual queda pendent
abans d’acceptar els resultats científicament.

## Validació del codi

Les proves amb dades sintètiques són a `tests/test_geography_mundlak.R`:

```r
source("tests/test_geography_mundlak.R", encoding = "UTF-8")
```

Comproven ponderació per població, reconstrucció de coordenades, recuperació
de pendents coneguts, derivades de les prediccions, contrasts amb covariància
no diagonal, etiquetatge within/between, signatura de les fonts, recuperació
d’exportacions i serialització del motor. Les proves de fórmules i transformació
de coordenades s’executen si hi ha `lme4` i `sf`.

Les proves no ajusten els models espacials grans ni accedeixen a les dades reals.
En la preparació s'han superat 107 comprovacions amb R 4.6.0 (webR),
lme4 2.0-6 i sf 1.1-3, incloses les comprovacions de fórmules i coordenades.
Cal validar la primera execució amb sdmTMB 1.1.0 a l’entorn Windows del projecte.
