# Després dels ambientals: extracció i demografia

Preparat el 20 d'octubre de 2026. Punt d'entrada: `post_environment_mundlak.R`.
Carregar el script només defineix funcions; no inicia cap procés.

## Què ha canviat respecte dels scripts anteriors

| Peça revisada | Dependència detectada / decisió |
| --- | --- |
| `env_mundlak_by_voltinism_qc_final_v3_resilient.R` | Quatre ajustos: onset i offset per uni/multivoltines; vuit components ambientals within/between. |
| `env_mundlak_wait_then_refine.R` | Els refinats es guarden a part. Cal consultar-los explícitament; els originals no són substituïts. |
| `env_mundlak_refine_v2_isolated.R` | La v2 substitueix la v1 per a l'offset univoltí i els dos onset; conserva l'offset multivoltí refinat de v1. Es comproven el controlador i els processos fills. |
| `extract_compare_offset_qc_plasticity_v1_3.R` | Extractor anterior d'offset: no admet la nova formulació ni l'onset separat. Es conserva com a historial. |
| `prepare_final_qc_abundance_inputs.R` | Preservava exactament l'onset antic. Aquesta regla no serveix per a la nova extracció. |
| `10f_mundlak_plasticity_demography_QC_FINAL_spatial_15km_RECOVERY.R` | Es reutilitzen la versemblança, optimització, checkpoints, estimació d'incertesa i contrasts. S'ha separat l'exportació perquè es pugui recuperar sense repetir l'ajust. |
| `geography_mundlak_by_voltinism_qc_final.R` | Ja prepara els quatre models geogràfics sobre fenologia anual observada. És una cua independent, no la font de plasticitats d'aquesta cadena. |

Els models demogràfics anteriors continuen disponibles i proporcionen els valors
inicials. Els nous predictors requereixen estimar de nou tots els paràmetres.
No es canvia la cua ambiental ni s'engega cap model durant aquesta preparació.

## Comprovar ara, sense carregar models

Des de l'R project de phenoIMPACT, amb el script actualitzat:

```r
source("post_environment_mundlak.R", encoding = "UTF-8")
pheno_post_env$preflight()
```

El preflight llegeix metadades i fitxers petits de comprovacions. Mostra l'estat
dels quatre ajustos, la versió de refinament seleccionada i si encara viu algun
worker ambiental, de v1, el controlador de v2 o un dels seus processos fills.
També detecta un fill viu encara que el controlador ja no existeixi. Si no es pot
comprovar un procés o les metadades de propietat són corruptes, atura la comprovació.
No carrega `fit_mundlak.rds`, els refinats ni els ajustos demogràfics.

La font ambiental predeterminada és:

```r
output/phenology_plasticity/environment_mundlak_by_voltinism_qc_final_15km/run_de0890c85d4a
```

La selecció és explícita:

- Sense cua v2, es conserva la selecció v1/original anterior.
- Amb `postfit_refinement_queue_v2_isolated`, l'offset univoltí i els dos onset
  esperen el resultat v2, encara que el seu procés no hagi començat. Un resultat
  v1 anterior no substitueix silenciosament una v2 pendent o fallida.
- L'offset multivoltí conserva el refinat de v1. Si la configuració v2 indica
  una carpeta `old_refinement` diferent, es respecta aquesta procedència.
- Cada refinat exigeix `fit_refined_final.rds`, `refinement_completed.rds`
  i comprovacions acceptades. Un resum de cua, per si sol, no valida un model.

Els fitxers seleccionats i `refinement_version` es conserven a
`environmental_sources.csv`. Abans de cada etapa es torna a comprovar que no hi
hagi processos ambientals/refinaments actius i que la selecció no hagi canviat.
Aquest comprovador no modifica ni atura les cues.

La condició ambiental manté el criteri existent: convergència 0, Hessiana positiva
i gradient màxim < 0,001. El valor 0,002998022 comunicat per a
`offset_multivoltine` **encara no passa**. Acabar la cua no implica passar aquesta porta.
El preflight no substitueix la revisió dels residus ambientals.

## Executar quan les fonts estiguin acceptades i el worker hagi acabat

```r
post_job <- pheno_post_env$run()
pheno_post_env$status()
```

La cadena fa:

1. Extracció dels quatre ajustos ambientals acceptats, un per procés.
2. Preparació de la demografia multivoltina, reajust i diagnòstics descriptius.
3. Preparació de la demografia univoltina, reajust i diagnòstics descriptius.
4. Combinació de les taules d'efectes i resum d'estat.

Abans de cada reajust, una avaluació del model demogràfic antic sobre les seves
dades antigues comprova la compatibilitat del motor TMB. No s'optimitza el model antic.
Cada etapa pesada corre en un procés separat per alliberar tota la memòria després.
L'ordinador ha de continuar encès i sense entrar en suspensió.

Opcionalment, preparar i inspeccionar la mostra abans d'ajustar:

```r
post_job <- pheno_post_env$run(mode = "prepare")
# Quan aquesta preparació hagi acabat:
post_job <- pheno_post_env$run(mode = "fit")
```

La segona crida reutilitza l'extracció i la preparació. Repetir `run()` després d'una
interrupció recupera les etapes completades i els checkpoints compatibles.
Un canvi de fonts, codi, versions de paquets o configuració de l'optimitzador crea
un directori nou; no reutilitza resultats amb predictors diferents.
Un ajust acabat amb problemes numèrics es conserva per revisar-lo, sense bucle de reajustos.

Si hi ha més d'una execució demogràfica anterior completa, el preflight mostra les
candidates i demana seleccionar-ne una explícitament amb `demography_run`:

```r
post_job <- pheno_post_env$run(
  demography_run = "E:/.../mundlak_geography_plasticity_qc_final_spatial_15km/run_..."
)
```

També admet `project_root`, `environment_run`, `minimum_free_ram_GiB` (20 per defecte),
`minimum_free_disk_GiB` (40) i `background = FALSE` per executar-lo en primer pla.
Els mínims de memòria i disc són comprovacions prèvies, no estimacions del pic de consum.
S'utilitzen les biblioteques R de la sessió que llança la cua, incloent sdmTMB 1.1.0.
Cal tenir `callr`, `ps`, `digest`, `TMB`, `Matrix` i, per als diagnòstics espacials, `spdep`.

## Contracte científic

- Plasticitat contextual = pendent fix de l'anomalia + vuit contribucions
  d'interacció ambiental within/between + desviació de pendent de l'espècie.
  Les finestres i unitats de l'anomalia es mantenen: 60 dies per onset i 90 per offset,
  amb pendent cru en dies/°C.
- L'onset anual dels offsets es manté constant en calcular la derivada. És una
  plasticitat d'offset condicionada a onset, no l'efecte tèrmic total sobre l'offset.
- Es conserva el signe cru. L'orientació biològica positiva és avançament per onset
  i offset univoltí, i retard per offset multivoltí.
- Cada variable s'estandarditza sobre **totes les poblacions extretes del seu ajust**,
  abans de l'encreuament amb abundància. Ara també l'onset és específic de cada grup.
  Es documenten mitjana, SD i univers de l'estandardització.
- La mostra de partida és l'abundància FINAL-QC del Mundlak demogràfic anterior.
  Es retenen tots els anys d'una població si té les dues noves plasticitats.
  No es recuperen automàticament poblacions excloses de l'abundància anterior.
  Les pèrdues i les poblacions extretes fora d'aquella mostra queden documentades.
- Els components within/between de plasticitat i la geografia es recalculen
  sobre les poblacions demogràfiques retingudes, amb igual pes per població.
  No es tornen a estandarditzar els components.
- Es mantenen Gamma/log, pendents i intercepts independents d'espècie i localitat,
  AR1 poblacional amb salts d'any reals, camp espacial anual IID i malla de 15 km.
  Els índexs de les variables latents es remapen amb les identitats guardades.
- No es fa un LRT entre els antics i els nous predictors. Les comprovacions
  d'imbricació anteriors només s'apliquen quan els predictors són els mateixos.
  Els contrasts Wald within/between del model nou es continuen exportant.
- Les SE del component fix de plasticitat i del pendent aleatori s'exporten
  separadament. **No són una SE total de plasticitat**: no es disposa de tota la
  covariància fix–aleatori. Igual que en el flux demogràfic anterior, la incertesa
  de la plasticitat no es propaga a les SE demogràfiques. Això continua sent una
  limitació/possible anàlisi de sensibilitat, no una qüestió resolta per aquest script.

## Outputs i estat final

Tot queda sota `output/population_trends/post_environment_mundlak_15km/run_<signatura>/`:

- `environmental_sources.csv`, `source_manifest.csv`, `code/`: procedència i còpia del codi.
- `extraction/<analisi>/population_plasticity.csv`: les quatre extraccions.
- `population_overlap_<grup>.csv`, `plasticity_scaling_<grup>.csv`: mostra i escala.
- `prepared_mundlak_<grup>.rds`: entrada amb predictors i índexs actualitzats.
- `<grup>__mundlak/`: ajustos, checkpoints, efectes, contrasts i covariàncies.
- `MUNDLAK_*.csv`: taules combinades compatibles amb els noms d'output existents.
- `<grup>__mundlak/residual_diagnostics/`: distribució Gamma PIT condicionada,
  residus per valors ajustats, correlacions temporals, innovacions latents AR1 i
  Moran anual bilateral (499 permutacions; veïnat de fins a 8 localitzacions).

Els diagnòstics condicionats són comprovacions descriptives, no tests de bondat
d'ajust calibrats amb simulacions i reajust. No s'hi apliquen p-valors IID de
uniformitat/dispersió. Les localitzacions coincidents s'agreguen per al Moran.

L'estat final numèric favorable és `FITTED_RESIDUAL_REVIEW_PENDING`, **no** una
declaració que tot el bloc científic ja estigui tancat. Si falla la convergència,
queda `NUMERICAL_REVIEW_REQUIRED`; si hi ha un error d'execució, el registre i els
checkpoints es conserven.

La cua geogràfica existent es pot executar separadament quan hi hagi recursos;
aquest script no l'engega ni canvia la seva formulació. Les figures definitives
s'actualitzaran després de revisar els nous resultats.

## Verificació del codi

```sh
Rscript tests/test_post_environment.R
```

Proves amb dades simulades: les quatre derivades i orientacions; onset additiu;
mapatge d'espècies; detecció de covariàncies/modes desfasats; filtre de poblacions;
escala abans de l'encreuament; components within/between; índexs de localitat/espècie;
AR1 amb anys absents; selecció original/refinat; portes numèriques; exports antics
i nous; recuperació dels exports; residus condicionats.

Aquestes proves i la lectura/sintaxi R s'han executat amb R via WebR. No s'han
executat els ajustos reals ni processos Windows en aquest entorn: els fitxers RDS
grans i el motor TMB compilat són a l'ordinador de treball. La comprovació definitiva
de compatibilitat nativa es fa amb el replay abans de començar els nous ajustos.
