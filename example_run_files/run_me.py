#!/usr/bin/env python3
import time
import inputR2S
from rmgcat_to_sella.main import WorkFlow

workflow = WorkFlow()
workflow.gen_job_files()
workflow.execute()
