from balsam.launcher.dag import BalsamJob


def restart():
    fails = []
    for job_status in ['FAILED', 'RUN_TIMEOUT']:
        fails.append(BalsamJob.objects.filter(state=job_status))

    for fail in fails:
        BalsamJob.batch_update_state(fail, 'RESTART_READY')
