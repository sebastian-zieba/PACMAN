from pathlib import Path

from astropy.io import ascii

from .options import OPTIONS


def nice_fit_par(fit_par_file: Path) -> None:
    fit_par = ascii.read(fit_par_file, Reader=ascii.CommentedHeader)
    with fit_par_file.open('w', encoding='ascii') as f:
        col_size = []
        for i in range(len(fit_par.colnames)):
            data = [str(ii) for ii in fit_par[fit_par.colnames[i]]]
            data_len = len(max(data, key=len))
            title_len = len(fit_par.colnames[i])
            col_size.append(max(title_len, data_len) + 2)

        file_header = []
        for i in range(len(fit_par.colnames)):
            file_header.append(fit_par.colnames[i])
            if i == 0:
                file_header.append(col_size[i]-1)
            else:
                file_header.append(col_size[i])

        file_rows = []
        for i in range(len(fit_par)):
            file_row = []
            for ii in range(len(fit_par.colnames)):
                file_row.append(fit_par[fit_par.colnames[ii]][i])
                file_row.append(col_size[ii])
            file_rows.append(file_row)

        print("#{: <{}} {: <{}} {: <{}} {: <{}} {: <{}} {: <{}} {: <{}} {: <{}} {: <{}} {: <{}} {: <{}} {: <{}}".format(*file_header), file=f)
        for row in file_rows:
            print("{: <{}} {: <{}} {: <{}} {: <{}} {: <{}} {: <{}} {: <{}} {: <{}} {: <{}} {: <{}} {: <{}} {: <{}}".format(*row), file=f)
