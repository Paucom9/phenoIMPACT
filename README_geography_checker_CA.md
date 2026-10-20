# Comprovador geogràfic: refinaments ambientals v1/v2

Versió 1.1.0, 20 d’octubre de 2026, de `geography_mundlak_by_voltinism_qc_final.R`.

El comprovador identifica la font ambiental acceptada abans de permetre els quatre models geogràfics. No modifica els ajustos ambientals, les cues de refinament ni els seus fitxers de resultats. Carregar el script i executar `preflight()` només llegeix metadades petites i consulta recursos i processos.

```r
source("geography_mundlak_by_voltinism_qc_final.R", encoding = "UTF-8")
pf <- pheno_geo_mundlak$preflight()
```

La taula `pf$source` mostra `kind`, `refinement_version`, `status`, `ready` i els camins seleccionats. `pf$sources_ready` indica si les quatre fonts són acceptables; `pf$ready` també exigeix que no hi hagi processos coneguts actius i que es compleixin els límits de recursos per defecte: 35 GiB de RAM disponible i 50 GiB de disc. `run()` torna a verificar les condicions amb els límits que rep.

## Selecció de fonts

- Sense refinaments, es comprova l’ajust original.
- Si hi ha v1 per a un model, té prioritat sobre l’original.
- Si existeix la cua v2, `offset_univoltine`, `onset_multivoltine` i `onset_univoltine` queden assignats a v2, encara que la seva carpeta individual no s’hagi creat.
- `offset_multivoltine` es conserva de v1. Es respecta `old_refinement` de la configuració v2 i es verifica que les configuracions pertanyin al mateix run ambiental.
- Un refinament pendent, incomplet o rebutjat no permet recuperar silenciosament un ajust anterior.

L’acceptació exigeix el marcador de finalització, el fitxer d’ajust, les comprovacions numèriques i l’entrada original preparada. Els estats admesos són `FIT_OK`, `FIT_OK_REFINED` i `FIT_OK_AFTER_REBUILD`, amb codi de convergència 0, Hessiana positiva, `numerical_checks_ok = TRUE` i gradient màxim absolut estrictament inferior a 0,001. Un resultat `OK` del procés aïllat no substitueix aquests controls.

Els models geogràfics continuen utilitzant les entrades originals preparades: observacions anuals, files, finestres climàtiques, estructura aleatòria i malla de 15 km. El canvi afecta la validació de la font; no canvia les fórmules ni extreu coeficients ambientals.

## Processos i canvis de font

Es comproven els controladors i els fitxers `worker_owner.rds`, `fit_owner.rds` i `isolated_worker_owner.rds`, inclòs un procés aïllat viu amb el controlador aturat i una font v1 configurada fora de les carpetes habituals. La identitat del procés exigeix PID i hora de creació, amb tolerància absoluta inferior a un mil·lisegon. Les metadades corruptes o els errors d’accés bloquegen la comprovació.

La selecció acceptada queda registrada a `environmental_source_selection.csv` en un futur run geogràfic. Abans de cada ajust es comprova que els ambientals continuïn inactius i que la selecció i les metadades dels fitxers no hagin canviat. No es carreguen els ajustos ambientals per fer aquestes comprovacions.

Per iniciar els models cal cridar explícitament `pheno_geo_mundlak$run()` quan les fonts siguin acceptables, els processos hagin acabat i hi hagi recursos. Aquest document no programa ni inicia cap ajust.

## Validació

```text
Rscript --vanilla tests/test_geography_checker.R --require-native
```

Els casos sintètics comproven originals, v1/v2, prioritat i procedència, fonts incompletes, llindars numèrics, processos aïllats, bloqueig abans d’escriure o iniciar treball, canvis de font i serialització del controlador. Les comprovacions natives utilitzen el procés de prova existent; no creen ni aturen processos de model.
