# Demografia: sensibilitat amb almenys 10 anys observats per població

El fitxer `demography_min10_observed_years.R` defineix `pheno_demo_min10`.
Carregar-lo amb `source()` no inicia cap procés. El mode predeterminat de
`run()` és `prepare`: els ajustos requereixen `mode = "fit"` explícit.

## Criteri i mostra

Població = `SPECIES × SITE_ID`. Es requereixen **almenys 10 valors diferents de
`YEAR` amb abundància vàlida en les entrades demogràfiques actuals**. Poden tenir
salts; 10 observacions consecutives són suficients encara que últim−primer = 9.
No es tracta del filtre inicial de transsectes amb 10 anys de visites.

| Grup | Observacions abans → després | Poblacions | Espècies | Transsectes |
| --- | ---: | ---: | ---: | ---: |
| Multivoltins | 152.971 → 117.796 | 15.643 → 7.475 | 73 → 64 | 1.792 → 1.694 |
| Univoltins | 304.056 → 251.015 | 27.026 → 15.146 | 163 → 129 | 1.787 → 1.678 |

## Contracte científic i dependències

La referència és `output/population_trends/post_environment_mundlak_15km/run_d984406cdf33`.
Es reutilitzen les seves estimacions puntuals de plasticitat i les còpies de codi
guardades dins de `code/`, amb comprovació d'identitat dels inputs i ajustos.
Aquestes dades locals i el DLL `trend_spde` són necessaris; GitHub conté el codi,
no els ajustos o les dades. No cal tornar a ajustar els models ambientals.

Es conserva Gamma/log, les pendents i intercepts d'espècie i transsecte, AR1 amb
salts reals, camp espacial anual IID i malla de 15 km. Es mantenen els anys del
camp i el centrat temporal originals. Es remapen espècies, transsectes i AR1;
les mitjanes entre espècies i les desviacions dins d'espècie es recalculen sobre
les poblacions retingudes, amb el mateix pes per població que al flux original.
**Els z de plasticitat i la seva SD de referència es conserven** per comparar
efectes en unitats comunes. Tots els paràmetres del model es tornen a estimar.
La propagació de la incertesa de la plasticitat queda fora per decisió de l'usuari.

Calen les biblioteques R de l'ajust de referència (`sdmTMB 1.1.0`, `TMB 1.9.17`,
`Matrix 1.7-3` tal com es registren al `run_config.rds`), a més de `digest`, `ps`,
`callr` i `spdep` per als diagnòstics espacials. El preflight compara les versions
registrades literalment amb les instal·lades.

## Preparar sense ajustar

Executa des de l'arrel del projecte en la sessió R habitual:

```r
source("demography_min10_observed_years.R", encoding = "UTF-8")
pf <- pheno_demo_min10$preflight()
job10 <- pheno_demo_min10$run(mode = "prepare")
pheno_demo_min10$status(job10)
```

La preparació pot fer-se mentre continua el geogràfic: només llegeix els ajustos
demogràfics petits i prepara entrades separades. No crea cap objectiu TMB ni
avalua el DLL. `preflight()` llegeix les entrades per comptar anys, però no inicia
processos ni escriu fitxers.

## Iniciar els dos ajustos quan hi hagi recursos

```r
pf <- pheno_demo_min10$preflight()
stopifnot(isTRUE(pf$ready_to_fit))
job10 <- pheno_demo_min10$run(mode = "fit")
pheno_demo_min10$status(job10)
```

Des d'un llançador que elimina tots els descendents quan acaba, mantén una
sessió de terminal oberta i executa `pheno_demo_min10$run(mode = "fit",
background = FALSE)`. Aquesta crida queda ocupada fins que acaba la cua; el
terminal persistent en conserva el procés. El mode de segon pla és per a
llançadors que permeten processos independents i desactiva explícitament
`cleanup`, `cleanup_tree` i `poll_connection`.

Un controlador executa un ajust alhora, en processos nous que alliberen la memòria
entre etapes. Detecta propietaris de processos de les altres cues del projecte i
refusa iniciar si n'hi ha de vius. Admissió predeterminada: 35 GiB RAM i 40 GiB de
disc disponibles. Es comproven recursos abans de cada ajust; **no hi ha una
protecció contínua que mati l'ajust si la memòria baixa després d'iniciar-lo**.
No es modifiquen ni s'aturen les altres cues.

La recuperació utilitza el motor provat existent: checkpoint de l'optimització
aproximadament cada 90 s quan hi ha una avaluació, paràmetres optimitzats desats
abans de la Hessiana i reutilització dels ajustos acabats. Tornar a executar
`run(mode="fit")` recupera la mateixa signatura si el codi, fonts i versions són
iguals. Un canvi de codi o fonts produeix un directori nou. Abans del primer ajust
de cada grup es verifica la versemblança antiga sobre les dades antigues.

## Outputs i revisió

Tot queda a `output/population_trends/demography_min10_observed_years_15km/run_<signatura>/`:

- `sample_before_after.csv`, `population_selection_<grup>.csv`: selecció exacta.
- `prepared_mundlak_<grup>.rds`: entrades filtrades i remapades.
- `<grup>__mundlak/`: ajustos, checkpoints, efectes i diagnòstics residuals.
- `MUNDLAK_*.csv`: efectes individuals, contrasts dins/entre, tests conjunts i
  balanç onset–offset, amb el format del flux de referència.
- `comparison_with_baseline.csv`: coeficients, IC95%, canvi de signe, efecte
  percentual i canvi d'error estàndard davant la referència.
- `review_fit_checks.csv`, `review_interpretation.txt`: comprovacions i límits.
- `comparison_residual_screens.csv`, `comparison_onset_offset_balance.csv`:
  diagnòstics residuals resumits i contrasts de balanç abans/després.
- `source_manifest.csv`, `run_config.rds`, `code/`: procedència i codi immutable.

La revisió s'exporta en acabar els dos grups; també es pot regenerar amb:

```r
pheno_demo_min10$review(job10)
# Després de reiniciar R:
# job10 <- list(out = "E:/.../demography_min10_observed_years_15km/run_<signatura>")
```

Cal comparar sobretot signe, magnitud i IC de l'onset univoltí dins d'espècie,
l'offset entre espècies i els contrasts d'allargament. Un canvi també pot reflectir
la pèrdua d'espècies i de cobertura geogràfica. No s'apliquen LRT entre mostres
diferents ni es classifiquen aquests ajustos per diferències d'AIC. `FIT_OK` és
acceptació numèrica; la revisió residual continua diferenciada.

## Verificació sense ajustar

```sh
Rscript --vanilla tests/test_demography_min10.R
Rscript --vanilla tests/test_demography_min10.R --real
```

El primer test usa dades sintètiques i comprova el criteri d'anys observats,
la pèrdua d'espècies/transsectes, els índexs no alfabètics, els salts AR1, la malla,
les escales i les signatures. El segon afegeix les entrades locals reals i verifica
els recomptes anteriors; no ajusta ni avalua la versemblança.
