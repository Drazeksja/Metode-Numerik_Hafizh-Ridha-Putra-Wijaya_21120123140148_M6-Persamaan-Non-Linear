# Jacobi_48.py
import math

def jacobi_48(x0, y0, tol=1e-6, max_iter=200):
    x, y = x0, y0
    print("=== Jacobi (NIM 48: g1A & g2A) ===")
    print(f"{'i':<4}{'x':<14}{'y':<14}{'Δx':<12}{'Δy':<12}")
    print(f"{0:<4}{x:<14.6f}{y:<14.6f}{0.0:<12.6f}{0.0:<12.6f}")
    for i in range(1, max_iter+1):
        # g1A and g2A using old x,y (Jacobi)
        if 10 - x*y < 0:
            print("sqrt domain negatif -> berhenti")
            return x, y
        x_new = math.sqrt(10 - x*y)         # g1A
        y_new = ((57 - 3 * x * (y**2)))     # <-- note: original g2A as rearrangement? 
        # BUT typical g2A used in problem statement: y_new = 57 - 3*x*y^2 (direct),
        # which is just algebraic rearrangement; we will use that form:
        y_new = 57 - 3 * x * (y**2)         # g2A

        dx, dy = abs(x_new - x), abs(y_new - y)
        print(f"{i:<4}{x_new:<14.6f}{y_new:<14.6f}{dx:<12.6f}{dy:<12.6f}")
        if dx < tol and dy < tol:
            print("Konvergen pada iterasi ke-", i)
            return x_new, y_new
        x, y = x_new, y_new
    print("Tidak konvergen dalam batas iterasi.")
    return x, y

# run example
if __name__ == "__main__":
    jacobi_48(1.5, 3.5)
