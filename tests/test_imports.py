import importlib
import pytest

MODULES = [
    'adapter',
    'BCTool_GUI',
    'Bivisu_algorithm',
    'Chen_Church_algorithm',
    'GO_assessment',
    'ISA_algorithm',
    'LAS_algorithm',
    'OPSM_algorithm',
]

DEPENDENCIES = {
    'BCTool_GUI': ['streamlit', 'pandas', 'openpyxl', 'matplotlib'],
    'GO_assessment': ['goatools'],
}

@pytest.mark.parametrize('module', MODULES)
def test_imports(module):
    for dep in DEPENDENCIES.get(module, []):
        pytest.importorskip(dep)
    importlib.import_module(module)
