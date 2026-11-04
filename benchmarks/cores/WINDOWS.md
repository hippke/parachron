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
Clang compiler (both are required: Parachron does not build with the Microsoft compiler alone):

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

## 5. Download and run the test

Paste these two commands (one at a time). The first downloads the test script into your user folder:

```powershell
cd $HOME; curl.exe -LO https://raw.githubusercontent.com/hippke/parachron/main/benchmarks/cores/run_windows.ps1
```

The second runs it:

```powershell
powershell -ExecutionPolicy Bypass -File .\run_windows.ps1
```

The script does everything else on its own, in a new folder `parachron-test` in your user folder, separate from anything
else on your PC:

1. it installs Parachron there (this compiles it and takes 2–5 minutes; warnings during compilation are normal),
2. it runs a quick test (about 2 minutes),
3. it runs the full test (5–10 minutes).

Please do not use the PC while the tests run, and let them finish. The fans may get loud; that is expected. The PC does
not go to sleep during the test. At the end the script prints `Finished.` and opens the folder with the results.

(`-ExecutionPolicy Bypass` allows this one script to run; it does not change any setting of your PC.)

## 6. Send us the result

The test writes files `parachron_bench_<name of your PC>_<date>.txt` into the folder `parachron-test` in your user
folder (the script opens it at the end; or type `explorer $HOME\parachron-test`). Please send us **all** of them (from
the quick and the full run). They contain your processor, the number
of cores, the Windows version, the power plan and the timings; no personal data except the computer's name, which you
can delete from the file and its name if you prefer.

## If something goes wrong

The script stops with a red `ERROR:` line that says what is missing:

- **`Python 3 was not found`**: install Python (step 2), close PowerShell, open it again, and repeat step 5.
- **`The installation of Parachron failed`**: most often "C++ Clang tools for Windows" was not ticked in step 3. Open the
  Visual Studio Installer, click **Modify**, tick it, then run step 5 again in a new PowerShell window. If that does
  not help, send us the file `setup_log.txt` from the `parachron-test` folder.
- **The test stops with an error**: send us the `parachron_bench_...txt` files anyway; the error is in them.

To run only the quick test (for example, to check that everything works):

```powershell
cd $HOME; powershell -ExecutionPolicy Bypass -File .\run_windows.ps1 -Quick
```

Running the script again later is fine: it reuses the installation and adds new result files.

To remove everything afterwards, delete the folder `parachron-test` and the file `run_windows.ps1` in your user folder.
Python and the Build Tools can be uninstalled under Settings → Apps.

What the test measures is described in [`README.md`](README.md#what-the-test-measures-for-us).
