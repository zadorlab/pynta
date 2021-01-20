#!/usr/bin/env python3
from pynta.main import WorkFlow


def run():
    # instantiate a WorkFlow() class
    workflow = WorkFlow()
    # create all input files
    workflow.gen_job_files()
    # execute the workflow
    workflow.execute_all()


if __name__ == '__main__':
    run()
