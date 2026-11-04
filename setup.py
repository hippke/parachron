# This file exists solely to configure the dynamic C extension.
# All project metadata lives in pyproject.toml.

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from glob import glob
import platform
import struct
import os
import subprocess
import sys
import sysconfig
import shutil


##### Windows: Parachron is compiled with clang-cl (MEP needs C11 atomics and C99 complex numbers, which the Microsoft
##### compiler lacks). It comes with Visual Studio ("C++ Clang tools for Windows") or with LLVM.
def find_clang_cl():
    if sys.platform != "win32":
        return None
    p = shutil.which("clang-cl")
    if p:
        return p
    candidates = []
    vswhere = os.path.join(os.environ.get("ProgramFiles(x86)", r"C:\Program Files (x86)"),
                           "Microsoft Visual Studio", "Installer", "vswhere.exe")
    if os.path.exists(vswhere):
        try:
            out = subprocess.check_output([vswhere, "-all", "-prerelease", "-products", "*", "-property", "installationPath"])
            for inst in out.decode(errors="replace").splitlines():
                if inst.strip():
                    candidates.append(os.path.join(inst.strip(), "VC", "Tools", "Llvm", "x64", "bin", "clang-cl.exe"))
        except Exception:
            pass
    candidates.append(os.path.join(os.environ.get("ProgramFiles", r"C:\Program Files"), "LLVM", "bin", "clang-cl.exe"))
    for c in candidates:
        if os.path.exists(c):
            return c
    return None

CLANG_CL_MISSING = """
Parachron needs the clang-cl compiler on Windows, and it was not found.
Install it with the Visual Studio Installer: Visual Studio Build Tools (or Visual Studio) -> Modify ->
workload "Desktop development with C++" -> also tick "C++ Clang tools for Windows" -> Modify.
Then run pip again in a new terminal.
"""

##### Compile with AVX512 if we are on x86 with 64bit.
def is_x86_64():
    if sys.platform == "darwin":
        return False
    if sys.platform == "win32":
        if not find_clang_cl():
            return False
    if sys.platform == "linux":
        import platform
        arch = platform.machine().lower()
        if not ("x86" in arch or "amd64" in arch):
            return False
        if "arm" in arch or "aarch" in arch:
            return False
    return struct.calcsize("P")*8 == 64

##### Git hash
try:
    ghash = subprocess.check_output(["git", "rev-parse", "HEAD"]).decode("ascii").strip()
    ghash_arg = f"-DGITHASH={ghash}"
except Exception:
    ghash_arg = "-DGITHASH=b06e3fcd9e50e2d9674f6f859ecd311091dd48b0" #GITHASHAUTOUPDATE

##### Link args
extra_link_args = []
if sys.platform == "darwin":
    cfg = sysconfig.get_config_vars()
    cfg["LDSHARED"] = cfg["LDSHARED"].replace("-bundle", "-shared")
    suffix = sysconfig.get_config_var("EXT_SUFFIX") or ".so"
    extra_link_args += [f"-Wl,-install_name,@rpath/libparachron{suffix}"]
if sys.platform.startswith("linux"):
    extra_link_args += ["-pthread", "-lrt"]   # mean-element parareal (src/mep.c): threads, semaphores, clock_gettime

##### Compile args
if sys.platform == "win32":
    extra_compile_args = [ghash_arg, "/wd4996", "/wd4244", "-DBUILDINGLIBREBOUND", "-D_GNU_SOURCE", "-DSERVER", "-D_CRT_NONSTDC_NO_WARNINGS=1",
                          "/clang:-ffp-contract=off", "-Wno-unused-function", "-Wno-unknown-pragmas"]
else:
    extra_compile_args = [ "-fstrict-aliasing", "-std=c99", "-Wno-unreachable-code", "-Wno-unknown-pragmas", ghash_arg, "-D_GNU_SOURCE", "-DSERVER", "-fPIC"]
    if os.environ.get("COVERAGE"):
        extra_compile_args += ["-O1", "-fprofile-arcs", "-ftest-coverage", "--coverage"]
        extra_link_args += ["-fprofile-arcs", "-ftest-coverage", "--coverage"]
    else:
        extra_compile_args.append("-O3")
    ##### Turn off floating point contractions for bitwise reproducibility
    if os.environ.get("FFP_CONTRACT_OFF"):
        extra_compile_args.append("-ffp-contract=off")
    ##### Optional shared-memory parallelisation (same as `make OPENMP=1` for the C library)
    if os.environ.get("REBOUND_OPENMP"):
        extra_compile_args += ["-fopenmp", "-DOPENMP"]
        extra_link_args += ["-fopenmp"]


class build_ext_avx512(build_ext):
    def build_extensions(self):
        os.makedirs(self.build_temp, exist_ok=True)
        if sys.platform == "darwin":
            if hasattr(self.compiler, 'linker_so'):
                for i, flag in enumerate(self.compiler.linker_so):
                    if flag == '-bundle':
                        self.compiler.linker_so[i] = '-shared'
        if sys.platform == "win32":
            # Use clang-cl instead of cl.exe with the Microsoft linker. The compiler is set after initialisation
            # (replacing the spawn function does not reach the compile calls of current setuptools).
            clang_cl = find_clang_cl()
            if clang_cl is None:
                raise RuntimeError(CLANG_CL_MISSING)
            if not getattr(self.compiler, "initialized", True):
                self.compiler.initialize()
            self.compiler.cc = clang_cl
            for opts in ("compile_options", "compile_options_debug"):
                if hasattr(self.compiler, opts):
                    setattr(self.compiler, opts, [o for o in getattr(self.compiler, opts) if o.upper() != "/GL"])
            print("parachron: compiling with " + clang_cl)
        if is_x86_64():
            if sys.platform == "win32":
                asm = os.path.join(self.build_temp, "integrator_whfast512_asm.obj")
                clang = os.path.join(os.path.dirname(find_clang_cl()), "clang.exe")
                cmd = [clang, "-c", "src/integrator_whfast512.s", "-o", asm]
                try:
                    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                except subprocess.CalledProcessError as e:
                    print(f"Assembly compilation failed: {e.stderr.decode()}", file=sys.stderr)
                    raise
            else:
                # Linux
                asm = os.path.join(self.build_temp, "integrator_whfast512_asm.o")
                self.compiler.spawn(["as", "--noexecstack", "-g", "-o", asm, "src/integrator_whfast512.s"])
            for ext in self.extensions:
                ext.extra_objects = (ext.extra_objects or []) + [asm]
        super().build_extensions()

##### C Extension
# The compiled C library (REBOUND C API, reb_* functions), named libparachron in the Python package so that it
# can be installed next to upstream REBOUND (librebound)
libreboundmodule = Extension(
    "libparachron",
    sources=sorted(glob("src/*.c")),
    include_dirs=["src"],
    depends=sorted(glob("src/*.s")),
    extra_link_args=extra_link_args,
    extra_compile_args=extra_compile_args,
)

setup(ext_modules=[libreboundmodule], cmdclass={"build_ext": build_ext_avx512})
