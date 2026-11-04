# Core-scaling benchmark (for testers, macOS)

This test measures how fast Parachron runs on your Mac with 1 core, then 2, 4, 8, ... up to all its cores. It takes
about 15 minutes of your time in total: about 5 minutes to install, and 5–10 minutes during which the test runs on
its own. At the end you send us one small text file. Nothing is sent anywhere automatically.

You type everything below into the **Terminal** app. Copy each grey box, paste it into Terminal (Cmd+V) and press
Return. Wait until the command has finished (the prompt, ending in `%` or `$`, appears again) before pasting the next
one.

## 1. Prepare the Mac

- Plug in the power adapter (on battery, macOS slows the processor down).
- Turn off Low Power Mode: System Settings → Battery → Low Power Mode → **Never**.
- Quit other apps, especially browsers, video calls and anything that is syncing or backing up.

## 2. Open Terminal

Press Cmd+Space, type `Terminal`, and press Return. A window with a text prompt opens.

## 3. Install Apple's developer tools

Parachron is compiled on your Mac during installation, so it needs Apple's free compiler:

```bash
xcode-select --install
```

A window appears: click **Install** and accept the license. This takes 5–15 minutes. If Terminal says
`command line tools are already installed`, that is fine: continue with step 4.

Check that Python works (it comes with these tools):

```bash
python3 --version
```

This should print `Python 3.9` or newer.

## 4. Install Parachron

These commands create a folder `parachron-test` in your home folder and install Parachron into it, separate from
anything else on your Mac:

```bash
mkdir -p ~/parachron-test && cd ~/parachron-test
```

```bash
python3 -m venv venv && source venv/bin/activate
```

The prompt now starts with `(venv)`.

```bash
pip install --upgrade pip
```

```bash
pip install parachron
```

The last command compiles Parachron and takes 1–3 minutes. It ends with `Successfully installed parachron-...`.
Warnings during compilation are normal.

## 5. Download the test

```bash
curl -LO https://raw.githubusercontent.com/hippke/parachron/main/benchmarks/cores/bench_cores.py
```

## 6. Run the test

First a quick check (about 2 minutes), to see that everything works:

```bash
python bench_cores.py --quick
```

If it ends with `Done. Please send the file ...`, run the full test:

```bash
caffeinate -i python bench_cores.py
```

(`caffeinate -i` keeps the Mac from going to sleep during the test.) Please do not use the Mac while it runs, keep a
laptop's lid open, and let it finish (5–10 minutes). The fans may get loud; that is expected.

## 7. Send us the result

The test writes a file `parachron_bench_<name of your Mac>_<date>.txt`. To find it, open the folder in Finder:

```bash
open ~/parachron-test
```

Please send us the files from **both** runs (the quick one and the full one). They contain the model of your Mac,
its processor, macOS version and the timings; no personal data except the computer's name, which you can delete from
the file and its name if you prefer.

## If something goes wrong

- **`xcrun: error: invalid active developer path`** or **`command not found: python3`**: step 3 has not finished. Run
  `xcode-select --install` again and wait until it is done.
- **`pip install parachron` fails**: save the full output and send it to us:

  ```bash
  pip install --no-cache-dir -v parachron > install_log.txt 2>&1
  ```

  The file `install_log.txt` is in `~/parachron-test`.
- **The test prints `WARNING: Python runs under Rosetta`**: your Python is an Intel version running in emulation, which
  would make the numbers meaningless. Tell us, and we will help.
- **The test stops with an error**: send us the `parachron_bench_...txt` file anyway; the error is in it.

To run the test again later (for example with the Mac on battery, for comparison), open Terminal and type:

```bash
cd ~/parachron-test && source venv/bin/activate && caffeinate -i python bench_cores.py
```

To remove everything afterwards, delete the folder `parachron-test` in your home folder (or run
`rm -rf ~/parachron-test`).

## What the test measures (for us)

`bench_cores.py` integrates the Solar System (Sun + 8 planets, WHFast512, dt = 6 d, GR, 17th-order corrector) and
reports:

1. **Single core**: ns per step of one simulation (median of 5 runs); directly comparable between machines.
2. **Throughput**: n independent simulations at the same time, one process each, for n = 1, 2, 4, ... physical cores
   (and the number of performance cores on CPUs with performance and efficiency cores). ns per step per simulation,
   and the efficiency relative to n × the single run.
3. **MEP**: one simulation on n cores with the default options (1600 slices of 60000 steps, 1.58 Myr), wall clock and
   speed-up against the serial time per step of section 1; the core counts are run up and then down, so that thermal
   throttling shows as a difference between the two passes.

The header records the CPU, the number of performance and efficiency cores, macOS version, compiler, power source,
Low Power Mode, and whether Python runs natively. It runs on Linux too (same output), for comparisons with Intel and
AMD machines.
