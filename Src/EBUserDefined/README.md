# EBUserDefined

User-replaceable embedded-boundary geometry for tools built with `USE_EB=TRUE`
and run with `eb2.geom_type = UserDefined` (`grad`, `curvature`,
`plotDisplacementSpeed`).

Replace `PeleLMeX_EBUserDefined.H` and `pelelmex_prob_parm.H` with the files of
your PeleLMeX case, rebuild (`make realclean` first) and add the geometry
parameters they read (typically `prob.*`) to the tool input file. The shipped
files are only an example geometry.

See the documentation page "Embedded boundaries (EB)"
(`Docs/source/basics/embeddedBoundaries.rst`) for details.
