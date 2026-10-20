# Sis sensibilitats filogenètiques Mundlak amb malla de 15 km

Versió 1.1.0 dels adaptadors i de la cua, 20 d’octubre de 2026.

| Família | Models | Component filogenètic |
| --- | --- | --- |
| Ambiental | onset_multivoltine, offset_multivoltine, onset_univoltine, offset_univoltine | Pendent d’anomalia climàtica per espècie |
| Demogràfica | univoltine, multivoltine | Pendent temporal per espècie |

Es conserven les files, els predictors, les escales, els efectes independents i la malla de cada font acceptada. Els ambientals mantenen Gaussian/ML, els vuit moderadors Mundlak i el control anual d’onset només als offsets; les finestres són de 60 dies per als dos onset i de 90 dies per als offsets. Els demogràfics mantenen Gamma/log, els 14 termes fixos, AR1 de població amb salts reals d’any, interceptes i pendents d’espècie i lloc i camp espacial anual.

La covariància conserva l’arrel original de l’arbre i diagonal unitària. Els predictors de plasticitat dels demogràfics actualitzats es mantenen fixos durant aquesta sensibilitat; la incertesa de la primera etapa no es propaga automàticament.

## Auditoria segura mentre corren els ambientals

```r
source("run_mundlak_phylogeny_weekend.R", encoding = "UTF-8")
audit <- audit_mundlak_phylogeny_weekend()
```

Aquesta crida només llegeix metadades petites i consulta processos. No carrega ajustos ni entrades preparades, no compila i no crea processos ni resultats. Mostra sis files i `ready = FALSE` si falten fonts, hi ha processos independents actius, alguna font no està acceptada o la procedència no concorda.

`start_mundlak_phylogeny_weekend()` també fa només aquesta auditoria per defecte (`check_only = TRUE`). En canvi, `run_environment_mundlak_phylogeny(check_only = TRUE)` i `run_demography_mundlak_phylogeny(check_only = TRUE)` són comprovacions de model que poden crear processos, carregar fitxers grans i construir exemples natius. No s’han d’utilitzar com a substitut de l’auditoria lleugera mentre hi ha ajustos en marxa.

## Fonts ambientals i refinaments

L’arrel per defecte és `output/phenology_plasticity/environment_mundlak_by_voltinism_qc_final_15km`, versió científica `environment_mundlak_by_voltinism_qc_final_15km_2.1_resilient`. L’onset agrupat anterior no és compatible amb els quatre models.

La selecció utilitza el comprovador geogràfic 1.1.0:

- V2 té prioritat per a `offset_univoltine`, `onset_multivoltine` i `onset_univoltine`, fins i tot abans que se’n creïn les carpetes individuals.
- `offset_multivoltine` es conserva de v1, respectant `old_refinement` de la configuració v2.
- Un refinament pendent, incomplet o rebutjat no recupera una versió anterior.
- Es comproven estat, identitat, fitxers i controls numèrics: convergència 0, Hessiana positiva i gradient màxim absolut inferior a 0,001.
- Els ajustos refinats es validen amb `source_signature` i `source_run`; les entrades preparades continuen sent les originals del mateix model i grup de voltinisme.

## Fonts demogràfiques actualitzades

L’arrel per defecte és `output/population_trends/post_environment_mundlak_15km`, amb versió `post_environment_mundlak_1.1`. La cua exigeix els dos ajustos completats i numèricament acceptats, no només una preparació o una carpeta existent.

Es verifiquen les quatre extraccions, la selecció ambiental v1/v2 desada, les metadades de les fonts, les signatures del run i els marcadors de preparació. La mida de la mostra prové de `prepared_<group>_completed.rds` i després es contrasta amb les dades i el contingut signat. No s’imposen les mides dels antics demogràfics, perquè la intersecció de poblacions amb plasticitat disponible pot haver canviat.

Es conserva el llindar numèric demogràfic `sqrt(g' V g) < 0,05`. La versemblança d’un demogràfic amb predictors o mostra actualitzats no es compara com a criteri d’acceptació amb la versemblança de la font anterior. Es mantenen les comprovacions internes del seu òptim, covariància, predictor i AIC.

L’auditoria individual permet consultar explícitament una font històrica 2.0.2. La cua de sis models exigeix fonts postambientals actualitzades i una mateixa procedència ambiental per als sis models. La revisió residual continua sent un pas separat; acceptació numèrica no equival a validació científica final.

## Execució posterior, explícita

Quan les sis fonts siguin acceptables, els processos independents hagin acabat i hi hagi recursos:

```r
job <- start_mundlak_phylogeny_weekend(
  check_only = FALSE,
  environment_source_run = "E:/.../environment_mundlak_by_voltinism_qc_final_15km/run_...",
  demography_source_run = "E:/.../post_environment_mundlak_15km/run_..."
)
watch_mundlak_phylogeny_weekend(job)
```

L’ordre és onset multivoltí → offset multivoltí → onset univoltí → offset univoltí → demogràfic univoltí → demogràfic multivoltí. Només hi ha un procés de model alhora. No s’ajusta automàticament cap font independent pendent. Cada llançament torna a comprovar processos i metadades, inclosos processos aïllats v2 amb un controlador aturat.

Es mantenen els mínims de RAM dels adaptadors: 30 GiB per als ambientals i 25 GiB per als demogràfics. La cua desa la selecció i les metadades de cada font; qualsevol canvi abans del següent model bloqueja aquell llançament.

## Fitxers i validació

Els adaptadors `11c_environment_mundlak_phylogeny_15km.R`, `10g_mundlak_demography_phylogeny_15km.R` i la cua utilitzen `mundlak_phylogeny_sources.R` i `geography_mundlak_by_voltinism_qc_final.R` 1.1.0. L’adaptador ambiental manté les dependències natives existents `pheno_phylogeny_qc_15km_v2.R` i `pheno_plasticity_phylogeny_pilot_15km.R`. No cal actualitzar paquets.

```text
Rscript --vanilla tests/test_mundlak_phylogeny_sources.R --require-native
```

Les proves utilitzen metadades sintètiques, dades petites en memòria i el procés de prova existent. No ajusten models ni compilen. La validació dels grans ajustos reals i les proves natives de versemblança queden per al llançament explícit posterior.
