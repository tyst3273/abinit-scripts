import Abinit_DFPT

workflow = Abinit_DFPT.workflow(input_file='INPUT',sbatch_job_file='job.sh')


workflow.write_abinit_phonon_input_files()

workflow.write_desktop_run_script()
#workflow.write_eagle_run_script()


