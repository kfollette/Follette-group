__author__ = 'JB'
import os

def kppPerDir(inputDir,obj_list,spec_path_list = None,outputDir = None, mute_error = True,compact_date_convention = None):

    if compact_date_convention is None:
        compact_date_convention = "GPIDATA"

    inputDir = os.path.abspath(inputDir)
    if compact_date_convention == "GPIDATA":
        compact_date=inputDir.split(os.path.sep)[-1].split("_")[0]

    # if outputDir is None:
    #     outputDir = inputDir
    # else:
    #     outputDir = os.path.abspath(outputDir)

    err_list = []
    for obj in obj_list:
        iterating = True
        while iterating:
            if not mute_error:
                iterating = obj.initialize(inputDir=inputDir,outputDir=outputDir,compact_date=compact_date)
                if obj.spectrum_iter_available() and spec_path_list is not None:
                    for spec_path in spec_path_list:
                        obj.init_new_spectrum(spec_path)
                        run(obj)
                else:
                    run(obj)
            else:
                try:
                    iterating = obj.initialize(inputDir=inputDir,outputDir=outputDir,compact_date=compact_date)

                    if obj.spectrum_iter_available() and spec_path_list is not None:
                        for spec_path in spec_path_list:
                            try:
                                obj.init_new_spectrum(spec_path)
                                run(obj)
                            except Exception as myErr:
                                err_list.append(myErr)
                                print("//!\\\\ "+obj.filename+"with spectrum "+spec_path+" in "+inputDir+" raised an Error.")
                    else:
                        try:
                            run(obj)
                        except Exception as myErr:
                            err_list.append(myErr)
                            print("//!\\\\ "+obj.filename+" in "+inputDir+" raised an Error.")
                except Exception as myErr:
                    err_list.append(myErr)
                    iterating = False
                    print("//!\\\\ "+obj.filename+" could NOT initialize in "+inputDir+". raised an Error.")


    return err_list


def run(obj):
    if not obj.check_existence():
        map = obj.calculate()
        obj.save()
    else:
        map = obj.load()

    return None