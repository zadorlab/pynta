import time
from typing import Any, List, Dict

from balsam.launcher.dag import BalsamJob
from balsam.core.models import ApplicationDefinition

from ase import Atoms
from ase.io import read, write
from ase.calculators.calculator import (
    Calculator, FileIOCalculator, all_changes
)


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
    command = '{exe} {args}'
    exe = None
    inpname = None
    outname = None
    args = None
    app = None

    preprocess = ''
    postprocess = ''
    description = ''

    ignore_fail = False

    def __init__(
        self,
        workflow: str,
        label: str = 'balsam',
        atoms: Atoms = None,
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
        self.job_kwargs = job_kwargs
        self.create_application()

    @classmethod
    def create_application(cls) -> None:
        """Creates a Balsam Application for this Calculator."""
        if cls.app is not None:
            return
        cls.app = ApplicationDefinition(
            name=cls.__name__,
            executable=cls.exe,
            preprocess=cls.preprocess,
            postprocess=cls.postprocess,
            description=cls.description,
        )
        cls.app.save()

    def write_input(
        self,
        atoms: Atoms,
        properties: List[str] = None,
        system_changes: List[str] = None
    ) -> None:
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        write(
            self.inpname.replace('PREFIX', self.label),
            atoms,
            **self.parameters
        )

    def calculate(
        self,
        atoms: Atoms = None,
        properties: List[str] = ['energy'],
        system_changes: List[str] = all_changes
    ) -> None:
        Calculator.calculate(self, atoms, properties, system_changes)
        job = BalsamJob(
            name=self.prefix,
            workflow=self.workflow,
            application=self.app.name,
            args=self.args.replace('PREFIX', self.prefix),
            **self.job_kwargs
        )

        self.directory = job.working_directory
        self.write_input(self.atoms, properties, system_changes)

        job.save()

        while True:
            job.refresh_from_db()
            if job.state in UNFINISHED_STATES:
                time.sleep(10)
                continue
            elif job.state in FAIL_STATES:
                if self.ignore_fail:
                    break
                else:
                    raise RuntimeError("Balsam job failed: " + job.state)
            elif job.state == 'JOB_FINISHED':
                break
            else:
                raise RuntimeError("Unknown Balsam job state: " + job.state)

        self.read_results()

    def read_results(self) -> None:
        out = read(self.outname.replace('PREFIX', self.label))
        self.out_calc = out.calc
        self.results = out.calc.results


class EspressoBalsamCalculator(BalsamCalculator):
    implemented_properties = ['energy', 'forces', 'stress', 'magmoms']
    exe = 'pw.x'
    inpname = 'PREFIX.pwi'
    outname = 'PREFIX.out'
    args = f'-in {inpname}'
    ignore_fail = True
