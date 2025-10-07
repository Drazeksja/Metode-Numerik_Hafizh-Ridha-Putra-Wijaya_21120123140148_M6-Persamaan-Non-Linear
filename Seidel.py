# Seidel_48.py
import math

def seidel_48(x0, y0, tol=1e-6, max_iter=200):
    x, y = x0, y0
    print("=== Seidel (NIM 48: g1A & g2A) ===")
    print(f"{'i':<4}{'x':<14}{'y':<14}{'Δx':<12}{'Δy':<12}")
    print(f"{0:<4}{x:<14.6f}{y:<14.6f}{0.0:<12.6f}{0.0:<12.6f}")
    for i in range(1, max_iter+1):
        # g1A using old x,y
        if 10 - x * y < 0:
            print("sqrt domain negatif -> berhenti")
            return x, y
        x_new = math.sqrt(10 - x * y)             # g1A
        # Use x_new when computing y (Seidel)
        # based on example in modul: y_{r+1} = ((57 - y_r) / (3 * x_{r+1}))**(1/3)
        denom = 3 * x_new
        if denom == 0:
            print("Denominator nol -> berhenti")
            return x_new, y
        inner = (57 - y) / denom
        # inner could be negative; cube root supports negative.
        y_new = math.copysign(abs(inner) ** (1/3), inner)

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
    seidel_48(1.5, 3.5)
