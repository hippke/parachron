# Core-scaling benchmark (for testers, Windows)

This test measures how fast Parachron runs on your PC with 1 core, then 2, 4, 8, ... up to all its cores. It takes
about 30 minutes of your time in total: about 20 minutes to install the tools (mostly waiting for downloads), and 5–10
minutes during which the test runs on its own. At the end you send us one small text file. Nothing is sent anywhere
automatically.

It needs Windows 10 or 11 on a 64-bit Intel or AMD processor (not ARM), and about 8 GB of free disk space for the
compiler.

## 1. Prepare the PC

- On a laptop, plug in the power adapter.
- Settings → System → Power (& battery): set **Power mode** to **Best performance**, and under **Screen and sleep**
  (or **Screen, sleep & hibernate timeouts**) set "when plugged in, put my device to sleep after" to **Never** for the
  duration of the test.
- Quit other programs, especially browsers, video calls, games, and anything that is syncing or updating.

## 2. Install Python

1. Download Python from https://www.python.org/downloads/ (the yellow "Download Python 3..." button).
2. Run the installer. On its first page, tick **"Add python.exe to PATH"**, then click **Install Now**.

If Python is already installed, you can skip this step.

## 3. Install the C compiler (Build Tools for Visual Studio)

Parachron is compiled on your PC during installation, so it needs Microsoft's free compiler tools together with the
Clang compiler:

1. Open https://visualstudio.microsoft.com/downloads/, scroll down to **Tools for Visual Studio**, and download
   **Build Tools for Visual Studio**.
2. Run the downloaded file. When the Visual Studio Installer shows the list of workloads, tick
   **Desktop development with C++**.
3. On the right side ("Installation details"), additionally tick **C++ Clang tools for Windows** (also called "C++
   Clang Compiler for Windows").
4. Click **Install** and wait until it has finished (5–20 minutes). Restart the PC if the installer asks you to.

If you already have Visual Studio, open the **Visual Studio Installer** from the Start menu, click **Modify** next to
your Visual Studio, and tick the same two items.

## 4. Open PowerShell

Press the Windows key, type `PowerShell`, and press Enter. A window with a text prompt opens.

You type everything below into this window. Copy each grey box, paste it into PowerShell (right-click or Ctrl+V) and
press Enter. Wait until the command has finished (the prompt, ending in `>`, appears again) before pasting the next
one.

## 5. Install Parachron

These commands create a folder `parachron-test` in your user folder and install Parachron into it, separate from
anything else on your PC:

```powershell
mkdir $HOME\parachron-test -Force; cd $HOME\parachron-test
```

```powershell
py -m venv venv
```

```powershell
.\venv\Scripts\python -m pip install --upgrade pip
```

```powershell
.\venv\Scripts\python -m pip install parachron
```

The last command compiles Parachron and takes 2–5 minutes. It ends with `Successfully installed parachron-...`.
Warnings during compilation are normal.

## 6. Download the test

```powershell
curl.exe -LO https://raw.githubusercontent.com/hippke/parachron/main/benchmarks/cores/bench_cores.py
```

## 7. Run the test

First a quick check (about 2 minutes), to see that everything works:

```powershell
.\venv\Scripts\python bench_cores.py --quick
```

If it ends with `Done. Please send the file ...`, run the full test:

```powershell
.\venv\Scripts\python bench_cores.py
```

Please do not use the PC while it runs, and let it finish (5–10 minutes). The fans may get loud; that is expected.

## 8. Send us the result

The test writes a file `parachron_bench_<name of your PC>_<date>.txt`. To find it, open the folder in Explorer:

```powershell
explorer $HOME\parachron-test
```

Please send us the files from **both** runs (the quick one and the full one). They contain your processor, the number
of cores, the Windows version, the power plan and the timings; no personal data except the computer's name, which you
can delete from the file and its name if you prefer.

## If something goes wrong

- **`py` is not recognized**: Python is not installed (step 2), or PowerShell was opened before the installation.
  Close PowerShell, open it again, and repeat from step 5.
- **`Parachron needs the clang-cl compiler on Windows, and it was not found`**: in step 3, "C++ Clang tools for Windows"
  was not ticked. Open the Visual Studio Installer, click **Modify**, tick it, and run the install command of step 5
  again in a new PowerShell window.
- **`Microsoft Visual C++ 14.0 or greater is required`** or **`link.exe` not found**: step 3 is missing or not finished.
- **The installation fails otherwise**: save the full output and send it to us:

  ```powershell
  .\venv\Scripts\python -m pip install --no-cache-dir -v parachron *> install_log.txt
  ```

  The file `install_log.txt` is in the `parachron-test` folder.
- **The test stops with an error**: send us the `parachron_bench_...txt` file anyway; the error is in it.

To run the test again later, open PowerShell and type:

```powershell
cd $HOME\parachron-test; .\venv\Scripts\python bench_cores.py
```

To remove everything afterwards, delete the folder `parachron-test` in your user folder. Python and the Build Tools can
be uninstalled under Settings → Apps.

What the test measures is described in [`README.md`](README.md#what-the-test-measures-for-us).
