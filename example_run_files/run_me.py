#!/usr/bin/env python3
from rmgcat_to_sella.main import WorkFlow

workflow = WorkFlow()
workflow.gen_job_files()
workflow.execute_all()
