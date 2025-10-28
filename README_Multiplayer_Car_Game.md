# Multiplayer Car Game (Python + Pygame)

> A small racing project where I built a **custom 2D physics engine** and a **UDP multiplayer layer** from scratch. Written as a second‑year CS student: clear, practical, and honest about what works and what I’d improve next.

**Project report (DOCX):** [51617_6024_Collin, Andrew (Multiplayer Car game)_project.docx](sandbox:/mnt/data/51617_6024_Collin, Andrew (Multiplayer Car game)_project.docx)

---

## What this project is
A side‑on racing game built with **Python** and **Pygame**. It features:
- **Deterministic physics** (frame‑rate independent integration).
- **Rigid body rotation** and **SAT collision detection** for rotated shapes.
- A basic **car model**: engine torque curve, wheel slip (Pacejka‑style), weight transfer, suspension via raycasts, and wheel rotation visuals.
- **Procedural terrain** (Perlin‑style noise) with a **parallax** background.
- A lightweight **UDP server** that relays player state and **lobbies** with ready‑up and start countdown.

This is a learning project. I focused on correctness and clarity first, then performance and UX.

---

## Quick start

### 1) Install prerequisites
```bash
# Python 3.10+ recommended
python -m venv .venv
# Windows
.venv\Scripts\activate
# macOS/Linux
source .venv/bin/activate

pip install pygame numpy noise
```

> If the repository includes a `requirements.txt`, you can install from that instead:
> ```bash
> pip install -r requirements.txt
> ```

### 2) Server (UDP relay)
Start the server first (by default it uses a single UDP port):
```bash
python server.py
```
> If you want a custom address/port, add simple flags (if present in your version), e.g.:
> ```bash
> python server.py --host 0.0.0.0 --port 9999
> ```
> Otherwise, set the values at the top of `server.py` and make sure clients match them.

### 3) Client
Run a client in another terminal (and repeat for additional players):
```bash
python main.py
```
- Enter a **username**.
- **Create** a lobby or **join** an existing one.
- When everyone is **Ready**, the server starts the countdown and the race begins.

### 4) Settings file
`main.py` expects a small `settings.txt` next to it. The project logic reads **line 2** for an assets path and **line 4** for the target FPS, each terminated by a semicolon. A minimal example that works with the default layout:

```
; Multiplayer Car Game settings
assets/ ; path to assets (images/, sound/)
; (unused placeholder line)
60 ; target FPS
```

If you see a prompt for the settings path on launch, provide the full path to your `settings.txt`.

---

## Controls (default)
- **D** – Gas / accelerate
- **A** – Brake / reverse
- **ESC** – Quit / back to menu

> Note: Controls are easy to change in code. Some experimental builds add gears and extra bindings.

---

## How multiplayer works (simple relay)
- **Protocol:** UDP
- **Server:** Receives client state (`x, y, rotation, wheel positions`) and relays it to players in the same lobby.
- **Client:** Packs and sends local car state at a steady tick; unpacks and renders remote cars.

Data messages are short tagged strings (e.g. `JS` for join server, `READY`, `RP` for “remote position”), which keeps the code readable while prototyping.

---

## Physics, rendering and UI
- **Integration:** Symplectic Euler; deterministic with fixed `dt`.
- **Collision:** Separating Axis Theorem (SAT) with simple response.
- **Car model:** Engine torque curve → drivetrain → wheel forces; slip ratio → grip; weight transfer; suspension via raycasts.
- **Terrain:** Chunked generation from noise; cheap to render; parallax layers for depth.
- **UI:** Fuel gauge, tachometer (RPM), timer, pedals indicator; lobby + join/create menus with ready‑up.

---

## Project structure (high‑level)
```
.
├── main.py          # client/game: menus, car, terrain, rendering, networking client
├── server.py        # UDP relay: lobbies, ready‑up, start, state fan‑out
├── images/          # sprites (car, wheels, backgrounds, UI)
├── sound/           # engine (planned) and SFX
├── settings.txt     # tiny config (see example above)
└── README.md        # this file
```

Folder names can vary by repository. The core idea is a single client script and a simple UDP server.

---

## Known limitations
- **Audio** is minimal; engine pitch mapping is a stub in some builds.
- **No car‑to‑car collisions** in networked races (visual only).
- **Single map**; terrain seeds are supported, but there’s no in‑game selector yet.
- **No persistence** (no accounts, coins, or leaderboards saved to disk).

---

## Roadmap (what I’d do next)
Short‑term
- Client‑side **interpolation/extrapolation** and basic lag compensation.
- **Authoritative server** toggle for anti‑cheat and fair starts.
- **Multiple maps** (seeded noise + background themes) and basic **car variants**.
- Improved **audio** (engine, finish, UI; simple mixer).

Medium‑term
- **Accounts + cloud save** (separate service), **leaderboards**, and daily seeds.
- **Replay/ghosts** recording for solo practice alongside real multiplayer.
- **Cross‑platform** packaging (Windows/macOS/Linux; explore mobile runtime options).

Long‑term
- **Tournament mode** and telemetry export.
- **Modding hooks** for cars and tracks.
- Optional **rollback netcode** for tighter races.

---

## Build notes
This is a teaching/portfolio project. I kept the code explicit (few abstractions, minimal dependencies) so it’s easy to read, profile, and extend.

If you try it and hit issues, check:
1) client/server **IP/port** match,
2) `settings.txt` path and FPS value,
3) your Python/pygame install.

---

## License
Personal/educational use. If you plan to reuse the physics/server code in your own project, add a license of your choice (MIT/BSD are fine) and credit the original author.

---

## Credits
Built by **Andrew Collin** as part of a self‑directed learning project in physics programming and networking. Thanks to the open articles and videos that helped me understand SAT, rotational dynamics, and basic UDP patterns.
