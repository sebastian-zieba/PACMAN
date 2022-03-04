import sys, os, time
sys.path.append('/home/zieba/Desktop/Projects/Open_source/PACMAN/')
import getopt
import pacman.reduction.s00_table as s00
import pacman.reduction.s01_horizons as s01
import pacman.reduction.s02_barycorr as s02
import pacman.reduction.s03_refspectra as s03
import pacman.reduction.s10_direct_images as s10
import pacman.reduction.s20_extract as s20
import pacman.reduction.s21_bin_spectroscopic_lc as s21
import pacman.reduction.s30_run as s30
from pacman.lib.update_meta import update_meta
from pacman.lib import sort_nicely as sn


def usage():
    cmd = sys.argv[0]
    sys.stderr.write('Usage: python %s OPTION\n\n' % os.path.basename(cmd))
    sys.stderr.write(
        'Allowed OPTION flags are:\n'
        '  --s00        reads in fits files and creates filelist.txt\n'
        '  --s01        downloads positions of HST during observations\n'
        '  --s02        corrects the MJD to BJD using the positions of HST\n'
        '  --s03        downloads the stellar spectrum and creates a reference spectrum with the bandpass of the grism\n'   
        '  --s10        determines the position of the source by looking at the direct image\n' 
        '  --s20        extracts the spectra\n' 
        '  --s21        bins light curves\n'
        '  --s30        fits models to the extracted light curve(s)\n'
        '  --workdir    sets the work directory (only important if user does not want to use the latest work dir)\n'
        '  --eventlabel name of the event; will be used for the naming of the new subdirectory (= work directory)\n'
        '  -h, --help   lists instructions for usage\n'
        '\n')
    sys.exit(1)


def main():
    #parses command line input
    try: opts, args = \
            getopt.getopt(sys.argv[1:],
                "hov", ["h", "help", "s00", "s01", "s02", "s03", "s10", "s20", "s21", "s30", "workdir=", "eventlabel="]
            )
    except getopt.GetoptError: usage()

    #defaults for command line flags
    run_s00 = False  # reads in fits files and creates filelist.txt
    run_s01 = False  # downloads positions of HST during observations
    run_s02 = False  # corrects the MJD to BJD using the positions of HST
    run_s03 = False  # downloads the stellar spectrum and creates a reference spectrum with the bandpass of the grism
    run_s10 = False  # determines the position of the source by looking at the direct image
    run_s20 = False  # extracts the spectra
    run_s21 = False  # bins light curves
    run_s30 = False  # fits models to the extracted light curve(s)

    print(opts)
    for o, a in opts:
        if o in ("-h", "--help"): usage()
        elif o == "--s00":
            run_s00 = True    
            print([i[0] for i in opts])
            print("--eventlabel" not in [i[0] for i in opts])
            if "--eventlabel" not in [i[0] for i in opts]:
                assert False, "please use --eventlabel='some_event_name' when running --s00"
            for oo, aa in opts:
                print(oo)
                if oo == "--eventlabel":
                    eventlabel = aa
                    print('eventlabel: ', eventlabel)
        elif o == "--s01": run_s01 = True
        elif o == "--s02": run_s02 = True
        elif o == "--s03": run_s03 = True
        elif o == "--s10": run_s10 = True
        elif o == "--s20": run_s20 = True
        elif o == "--s21": run_s21 = True
        elif o == "--s30": run_s30 = True
        elif o == "--workdir": workdir = a
        elif o == "--eventlabel": continue
        else: assert False, "unhandled option. Please seek --help"

    if run_s00:
        meta = s00.run00(eventlabel)
        workdir = meta.workdir + '/'

    if not run_s00:
        # list subdirectories in the run directory
        dirs = [f.path for f in os.scandir('.') if f.is_dir()]
        # saves times when these subdirectories were created
        dirs_times = [i[6:25] for i in dirs]
        # sort the times
        times_sorted = sn.sort_nicely(dirs_times)
        # most recent time
        recent_time = times_sorted[-1]
        # find the directory with that most recent time
        idx = 0
        for i in range(len(dirs)):
            if dirs[i][6:25] == recent_time:
                idx = i
        workdir = dirs[idx]
        #save the eventlabel which is in the directory name too
        eventlabel = workdir[26:]
        workdir = workdir + '/'
        print('workdir: ', workdir)
        print('eventlabel: ', eventlabel)

    if run_s01:
        update_meta(eventlabel, workdir)
        meta = s01.run01(eventlabel, workdir)

    if run_s02:
        update_meta(eventlabel, workdir)
        meta = s02.run02(eventlabel, workdir)

    if run_s03:
        update_meta(eventlabel, workdir)
        meta = s03.run03(eventlabel, workdir)

    if run_s10:
        update_meta(eventlabel, workdir)
        meta = s10.run10(eventlabel, workdir)

    if run_s20:
        update_meta(eventlabel, workdir)
        meta = s20.run20(eventlabel, workdir)

    if run_s21:
        update_meta(eventlabel, workdir)
        meta = s21.run21(eventlabel, workdir)

    if run_s30:
        update_meta(eventlabel, workdir)
        meta = s30.run30(eventlabel, workdir)


if __name__ == '__main__':
    main()
