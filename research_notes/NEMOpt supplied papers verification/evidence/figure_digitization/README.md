# Figure-digitization scripts (evidence)

These scripts produced the digitized benchmark values in the verification notes. The source PDFs are **not** in the repository (copyright); obtain them separately and pass the PDF path as the first argument.

- `vector_michel2007/` (PyMuPDF): curves extracted from the vector figures of Michel, Lopez-Pamies, Ponte Castañeda & Triantafyllidis (2007), JMPS 55(5):900–938, doi:10.1016/j.jmps.2006.11.006, read from the open HAL deposit hal-00111501 (submitted version). Stated accuracy about ±0.002 in log strain. See `../../michel_et_al_2007.md`.
- `raster_triantafyllidis2006/`: values digitized from the figures of Triantafyllidis, Nestorović & Schraad (2006), J. Appl. Mech. 73(3):505–515, doi:10.1115/1.2126695 (PDF supplied by the authors of NEMOpt). Uncertainties are given per value in `../../triantafyllidis_nestorovic_schraad_2006.md`.

The split between the two folders follows the method each track reported (vector extraction for Michel et al.; raster digitization for Triantafyllidis et al.). The scripts were written by the research tracks and were not re-run by the coordinator, because the source PDFs are not stored with the repository. Check digitized values against the journal versions before fixing benchmark tolerances.
