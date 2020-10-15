#!/usr/bin/env python3
from pynta.main import WorkFlow

workflow = WorkFlow()
workflow.gen_job_files()
workflow.execute_all()
