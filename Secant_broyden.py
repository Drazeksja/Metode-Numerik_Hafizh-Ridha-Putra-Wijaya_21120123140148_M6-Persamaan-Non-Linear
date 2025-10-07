# Broyden_48.py
import numpy as np

def F(xy):
    x, y = xy
    return np.array([x**2 + x*y - 10, y + 3*x*(y**2) - 57], dtype=float)

def broyden_48(x0, y0, tol=1e-6, max_iter=100):
    xy = np.array([x0, y0], dtype=float)
    x, y = xy
    # initial Jacobian approximate = analytical Jacobian at initial guess
    J = np.array([[2*x + y, x],
                  [3*(y**2), 1 + 6*x*y]], dtype=float)
    F_val = F(xy)
    print("=== Broyden (Secant multivar) ===")
    print(f"{'i':<4}{'x':<14}{'y':<14}{'Δx':<12}{'Δy':<12}")
    print(f"{0:<4}{x:<14.6f}{y:<14.6f}{0.0:<12.6f}{0.0:<12.6f}")
    for i in range(1, max_iter+1):
        try:
            s = np.linalg.solve(-J, F_val)  # step
        except np.linalg.LinAlgError:
            print("Jacobian singular")
            return xy
        xy_new = xy + s
        dx, dy = abs(s[0]), abs(s[1])
        print(f"{i:<4}{xy_new[0]:<14.6f}{xy_new[1]:<14.6f}{dx:<12.6f}{dy:<12.6f}")
        if max(dx, dy) < tol:
            print("Konvergen pada iterasi ke-", i)
            return xy_new
        F_new = F(xy_new)
        # rank-1 Broyden update
        J = J + np.outer((F_new - F_val - J @ s), s) / (s @ s)
        xy, F_val = xy_new, F_new
    print("Tidak konvergen.")
    return xy

# run example
if __name__ == "__main__":
    broyden_48(1.5, 3.5)
