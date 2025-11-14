function calcAcceleration() {
  const F = parseFloat(document.getElementById("force").value);
  const m = parseFloat(document.getElementById("mass").value);
  const r = document.getElementById("physResult");

  if (!m || m === 0) {
    r.textContent = "Mass must be nonzero.";
    return;
  }

  const a = F / m;
  r.textContent = `Acceleration: ${a.toFixed(2)} m/s^2`;
}
