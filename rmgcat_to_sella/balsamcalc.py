import time
from typing import Any, List, Dict
import socket
import atexit
from balsam.launcher.dag import BalsamJob, kill
from balsam.core.models import ApplicationDefinition

from ase import Atoms
from ase.io import read, write
from ase.calculators.calculator import (
    Calculator, FileIOCalculator, all_changes
)
from ase.calculators.socketio import SocketIOCalculator


UNFINISHED_STATES = [
    'CREATED',
    'AWAITING_PARENTS',
    'READY',
    'STAGED_IN',
    'PREPROCESSED',
    'RUNNING',
    'RUN_DONE',
    'POSTPROCESSED'
]


FAIL_STATES = [
    'RUN_TIMEOUT',
    'RUN_ERROR',
    'RESTART_READY',
    'FAILED',
    'USER_KILLED'
]


class BalsamCalculator(FileIOCalculator):
    # Command executed by Balsam
    exe = None

    # Naming scheme for the input file written by the calculator
    inp_name = None
    # ASE IO format of the input file
    inp_format = None
    # Naming scheme for the output file read by the calculator
    out_name = None
    # ASE IO format of the output file
    out_format = None

    # Extra calculation-specific arguments to provide to Balsam
    args = None

    # Extra information for the Balsam App
    preprocess = ''
    postprocess = ''
    description = ''

    # Ignore when Balsam jobs fail (needed for QE in socket-mode)
    ignore_fail = False

    # Balsam job object
    job = None

    def __init__(
        self,
        workflow: str,
        label: str = 'balsam',
        atoms: Atoms = None,
        job_args: str = None,
        job_kwargs: Dict[str, Any] = dict(),
        **kwargs: Any
    ) -> None:
        FileIOCalculator.__init__(
            self,
            restart=None,
            ignore_bad_restart_file=False,
            label=label,
            atoms=atoms,
            **kwargs
        )
        self.workflow = workflow
        self.job_args = job_args
        self.job_kwargs = job_kwargs
        atexit.register(self.kill_job)
        self.create_application()

    @classmethod
    def create_application(cls) -> None:
        try:
            app = ApplicationDefinition.objects.get(
                name=cls.__name__
            )
        except ApplicationDefinition.DoesNotExist:
            app, created = ApplicationDefinition.objects.get_or_create(
                name=cls.__name__,
                executable=cls.exe,
                preprocess=cls.preprocess,
                postprocess=cls.postprocess,
                description=cls.description,
            )
            if created:
                app.save()

    def format_args(self) -> str:
        args = self.args.replace('PREFIX', self.prefix)
        if self.job_args is not None:
            args = self.job_args + ' ' + args
        return args

    def write_input(
        self,
        atoms: Atoms,
        properties: List[str] = None,
        system_changes: List[str] = None
    ) -> None:
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        write(
            self.inp_name.replace('PREFIX', self.label),
            atoms,
            format=self.inp_format,
            **self.parameters
        )

    def create_job(self) -> BalsamJob:
        return BalsamJob(
            name=self.prefix,
            workflow=self.workflow,
            application=self.__class__.__name__,
            args=self.format_args(),
            **self.job_kwargs
        )

    def kill_job(self) -> None:
        if self.job is None:
            return
        if self.job.state in UNFINISHED_STATES:
            kill(self.job)

    def job_running(self, state: str) -> bool:
        if state in UNFINISHED_STATES:
            return True
        if state in FAIL_STATES:
            if self.ignore_fail:
                return False
            raise RuntimeError("Balsam job failed: " + state)
        if state == 'JOB_FINISHED':
            return False
        raise RuntimeError("Unknown Balsam job state: " + state)

    def calculate(
        self,
        atoms: Atoms = None,
        properties: List[str] = ['energy'],
        system_changes: List[str] = all_changes
    ) -> None:
        Calculator.calculate(self, atoms, properties, system_changes)
        self.job = self.create_job()
        self.directory = self.job.working_directory
        self.write_input(self.atoms, properties, system_changes)
        self.job.save()

        while self.job_running(self.job.state):
            time.sleep(10)
            self.job.refresh_from_db()

        self.read_results()

    def read_results(self) -> None:
        out = read(
            self.out_name.replace('PREFIX', self.label), format=self.out_format
        )
        self.out_calc = out.calc
        self.results = out.calc.results


class EspressoBalsam(BalsamCalculator):
    implemented_properties = ['energy', 'forces', 'stress', 'magmoms']
    exe = 'pw.x'
    inp_name = 'PREFIX.pwi'
    inp_format = 'espresso-in'
    out_name = 'PREFIX.out'
    out_format = 'espresso-out'
    args = f'-in {inp_name}'
    ignore_fail = True


class BalsamSocketIOCalculator(BalsamCalculator, SocketIOCalculator):
    job = None

    def __init__(
        self,
        workflow: str,
        label: str = 'balsam',
        atoms: Atoms = None,
        job_args: str = None,
        job_kwargs: Dict[str, Any] = dict(),
        **kwargs: Any
    ) -> None:
        SocketIOCalculator.__init__(self, port=0)
        self._port = self.server.serversocket.getsockname()[1]
        BalsamCalculator.__init__(
            self,
            workflow=workflow,
            label=label,
            atoms=atoms,
            job_args=job_args,
            job_kwargs=job_kwargs,
            **kwargs
        )

    def calculate(
        self,
        atoms: Atoms = None,
        properties: List[str] = ['energy'],
        system_changes: List[str] = all_changes
    ) -> None:
        Calculator.calculate(self, atoms, properties, system_changes)
        if self.job is None or not self.job_running(self.job.state):
            self.job = self.create_job()
            self.directory = self.job.working_directory
            self.write_input(self.atoms, properties, system_changes)
            self.job.save()
        SocketIOCalculator.calculate(self, atoms, properties, system_changes)


class EspressoBalsamSocketIO(BalsamSocketIOCalculator):
    exe = 'pw.x'
    inp_name = 'PREFIX.pwi'
    inp_format = 'espresso-in'
    out_name = 'PREFIX.out'
    out_format = 'espresso-out'
    args = f'--ipi HOSTNAME:PORT -in {inp_name}'
    ignore_fail = True

    def format_args(self) -> str:
        return (
            BalsamCalculator.format_args(self)
            .replace('HOSTNAME', socket.gethostname())
            .replace('PORT', str(self._port))
        )
