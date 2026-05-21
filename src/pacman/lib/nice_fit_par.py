from pathlib import Path
from astropy.io import ascii


def nice_fit_par(fit_par_file: Path) -> None:
    fit_par = ascii.read(
        fit_par_file,
        format="commented_header",
        delimiter=r"\s",
        guess=False,
        fast_reader=False,
        fill_values=[("", "0")],
    )

    col_widths = []
    for col in fit_par.colnames:
        values = [str(value) for value in fit_par[col]]
        col_widths.append(max(len(col), max(len(v) for v in values)) + 2)

    with fit_par_file.open("w", encoding="ascii") as f:
        header = "#"
        header += "".join(
            f"{col:<{width}}"
            for col, width in zip(fit_par.colnames, col_widths)
        ).rstrip()
        print(header, file=f)

        for row in fit_par:
            line = "".join(
                f"{str(row[col]):<{width}}"
                for col, width in zip(fit_par.colnames, col_widths)
            ).rstrip()
            print(line, file=f)
            