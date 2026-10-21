# Rescat separat de l'offset multivoltí

Versió 1.0.0, 20 d'octubre de 2026. Script: `env_mundlak_rescue_offset_multivoltine.R`.

Parteix de `fit_refined_final.rds` de v1 per a `offset_multivoltine`. Respecta `old_refinement` de la configuració v2. Conserva dades, fórmula, malla de 15 km, grup de voltinisme, control d'onset, efectes aleatoris, Gaussian/ML, mapes i límits dels paràmetres. No modifica v1, v2, els scripts ambientals ni la selecció de fonts postambientals.

## Ús segur ara

```r
source("env_mundlak_rescue_offset_multivoltine.R", encoding = "UTF-8")
audit <- audit_offset_multivoltine_rescue()
# També: run_offset_multivoltine_rescue() — el mode per defecte és "audit".
```

L'auditoria només llegeix metadades petites, versions dels paquets, espai disponible i processos. No obre ajustos ni dades preparades, no compila i no inicia optimitzadors ni processos fills. `source()` només defineix funcions.

Els modes pesants es bloquegen abans de carregar el model si hi ha qualsevol procés independent detectat, inclosos el controlador v2 i fills aïllats sense controlador. No hi ha cap opció per saltar aquesta protecció. Un error llegint propietaris o comprovant processos també bloqueja. Requereix el comprovador geogràfic 1.1.0 del repositori.

## Execució futura explícita

Només quan no hi hagi processos independents actius i els recursos siguin suficients, en una consola PowerShell:

```powershell
# Diagnòstic pesant: reconstrueix l'objectiu i avalua gradients, sense nlminb.
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' --vanilla `
  'E:\phenoIMPACT project\code\phenoIMPACT\env_mundlak_rescue_offset_multivoltine.R' `
  --mode=diagnose --execute

# Rescat: pot executar nlminb i calcular la covariància final.
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' --vanilla `
  'E:\phenoIMPACT project\code\phenoIMPACT\env_mundlak_rescue_offset_multivoltine.R' `
  --mode=rescue --execute
```

Sense `--execute`, qualsevol mode pesant es rebutja. Sense arguments, Rscript fa només l'auditoria. Les funcions també exigeixen una sessió no interactiva per als modes pesants; no s'han d'executar dins la sessió RStudio de v2.

Es poden especificar `--project-root=...` i `--source-run=...`. Cada intent crea una carpeta nova dins `output/phenology_plasticity/offset_multivoltine_rescue/`. Un bloqueig `rescue.lock` impedeix dos rescats simultanis; si queda després d'una interrupció abrupta, s'ha d'inspeccionar el PID, la data de creació i la carpeta de resultats abans de retirar-lo.

El diagnòstic pesant pot consumir molta RAM i temps encara que no optimitzi: ha de llegir el gran ajust v1 i reconstruir el likelihood. No s'ha executat aquesta part sobre el model real durant la preparació del script.

## Diagnòstics i estratègia

- Desa el missatge, els controls, els comptadors i el gradient/covariància originals disponibles.
- Reconstrueix amb `do_fit = FALSE` i exigeix concordança de l'objectiu amb v1 (tolerància absoluta 0,00001).
- Avalua objectiu i gradient tres vegades al mateix vector. Exporta gradients amb signe, absoluts, noms, índexs únics, termes dels coeficients fixos, límits, distància als límits i variació entre repeticions.
- Si objectiu o gradient no són reproduïbles dins les toleràncies declarades, s'atura amb `UNSTABLE_GRADIENT`.
- Utilitza `nlminb` directament amb `rel.tol = 1e-12`, `x.tol = 1e-10`, `iter.max = 100` i `eval.max = 200`. Són controls explícits d'una primera prova; no garanteixen convergència. No activa Newton, no canvia l'escalat i no altera el llindar d'acceptació.
- Registra cada avaluació d'objectiu/gradient, el temps i els comptadors de cada passada. No recalcula `sdreport()` a cada passada.
- Per defecte hi ha com a màxim dues passades, i s'atura després de la primera sense progrés mesurable.

Es considera progrés si baixa l'objectiu més que el màxim de la tolerància numèrica, `1e-7` i `1e-12 * max(1, abs(objectiu))`, **o** si el gradient màxim baixa més que el màxim de `1e-5` i l'1% del gradient anterior. Una baixada clara del gradient pot justificar continuar encara que l'objectiu gairebé no canviï. La repetició exacta de les passades 3/4 de v1 s'atura amb `STAGNATION`.

Un objectiu pitjor més enllà del marge de soroll es rebutja. Un gradient pitjor sense progrés en l'objectiu tampoc es conserva. `retained` només significa que es conserva el candidat; no significa acceptació numèrica.

## Límits i aturada

Valors per defecte: 40 GiB de RAM disponible, 40 GiB de disc, 1 hora per passada d'optimització i 6 hores totals. Són límits cooperatius: es comproven entre avaluacions i operacions. Una crida nativa, una lectura gran o el càlcul de `sdreport()` que ja ha començat pot superar el temps abans de retornar. No hi ha cap procés que mati R automàticament.

Crear un fitxer buit `STOP_REQUESTED` dins la carpeta d'aquest intent demana una aturada cooperativa del rescat. No té cap efecte sobre v2. Els checkpoints ja escrits es conserven.

Per modificar controls futurs, es pot desar una llista RDS i passar `--control-file=...`. Exemple de preparació, que no executa res:

```r
saveRDS(list(max_passes = 1L, maximum_pass_seconds = 1800),
        "offset_rescue_control.rds")
```

Les opcions desconegudes, no finites o invàlides es rebutgen. El llindar científic del gradient no és configurable des d'aquest script.

## Resultats i acceptació

Els fitxers petits inclouen `source_optimizer_and_covariance.rds`, diagnòstics inicials i per passada, `rescue_history.csv`, traces d'avaluacions, registres d'errors i checkpoints amb paràmetres conjunts. El checkpoint d'avaluació és l'últim punt avaluat, que pot no ser el millor; queda etiquetat com a no acceptat. Els checkpoints serveixen per inspecció i recuperació explícita futura; el script no els reprèn automàticament.

Només quan el gradient estable baixa de 0,001 amb convergència 0 es calcula la covariància final, amb `par.fixed` explícit i sense la gran covariància de tots els informes. Després es tornen a comprovar:

1. Convergència 0, gradient absolut màxim < 0,001 i repetibilitat.
2. Hessiana positiva i covariància finita positiva definida.
3. Coincidència del punt i del gradient entre `sdreport` i l'avaluació directa.
4. Objectiu no pitjor que el punt inicial més enllà del soroll permès.
5. Procedència i fitxers font inalterats, i absència d'altres processos actius.

El gradient escalat amb la covariància es desa com a diagnòstic complementari. No substitueix el criteri de gradient absolut.

Si passa, escriu `fit_rescued_final.rds`, `fit_checks_rescued.csv` i un resultat `FIT_OK_RESCUED_REVIEW_REQUIRED`. El fit conserva les dades i els paràmetres portables, però no desa l'objectiu natiu viu. **No queda registrat automàticament com a font** dels fluxos postambiental, geogràfic o filogenètic: aquesta integració i la revisió científica es faran explícitament després de validar el resultat.

Una aturada, error, estancament o Hessiana fallida mai produeix un resultat acceptable ni canvia marcadors v1/v2.

## Validació

```text
Rscript --vanilla tests/test_offset_multivoltine_rescue.R --require-native
```

Les proves utilitzen metadades sintètiques, objectius matemàtics de dues variables i un optimitzador simulat. Verifiquen estancament, progrés, límits, errors, checkpoints, gradients inestables, Hessianes invàlides i bloqueig per processos v2 aïllats. No executen `nlminb`, cap likelihood TMB, cap ajust real, cap compilació ni cap procés de model.
