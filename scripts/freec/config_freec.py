import importlib
import subprocess
import sys

def install_and_import(package_name, import_name=None, upgrade=False):
    if import_name is None:
        import_name = package_name
    
    try:
        module = importlib.import_module(import_name)
        print(f"'{import_name}' installed, version: {getattr(module, '__version__', 'Unknown')}")
        return module
    except ImportError:
        print(f"Not found '{import_name}', Installing '{package_name}'...")
        pip_args = [sys.executable, "-m", "pip", "install", package_name]
        if upgrade:
            pip_args.append("--upgrade")
        
        try:
            subprocess.check_call(pip_args)
            print(f"'{package_name}' Installed.")
            return importlib.import_module(import_name)
        except subprocess.CalledProcessError as e:
            print(f"Install '{package_name}' failed!: {e}")
            sys.exit(1)

import configparser
# from snakemake.script import snakemake
config = configparser.ConfigParser()
config.optionxform = str ## don not covert to lower case!!!


## config temple
config_tem = snakemake.input['ini_template']
## input BAM file,for WES, a target panel bed
tumor_bam = snakemake.input['Tum']
normal_bam = '' if "NC" not in snakemake.input.keys() else snakemake.input['NC']
## freec params
chrLenFile = snakemake.params['chrLenFile']
outputDir = snakemake.params['outputDir']
chrFiles = snakemake.params['chrFiles']
bed = snakemake.params['bed']
sambamba = snakemake.params['sambamba']
threads = str(snakemake.threads)
### output config for each sample
out_config = snakemake.output['config']

config.read(config_tem)

config['general']['chrLenFile'] = chrLenFile
config['general']['outputDir'] = outputDir
config['general']['chrFiles'] = chrFiles
config['general']['sambamba'] = sambamba
config['general']['SambambaThreads'] = threads
config['general']['maxThreads'] = snakemake.params['maxThreads']


# config['BAF']['makePileup'] = snakemake.params['snp_file']
# config['BAF']['fastaFile'] = snakemake.params['ref']
# config['BAF']['SNPfile'] = snakemake.params['snp_file']

config['sample']['mateFile'] = tumor_bam
config['control']['mateFile'] = normal_bam
if "NC" not in snakemake.input.keys():
    config.remove_section('control')


config['target']['captureRegions'] = bed

with open(out_config, 'w') as configfile:
    config.write(configfile)
