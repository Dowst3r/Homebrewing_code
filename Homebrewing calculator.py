import tkinter as tk
import tkinter.messagebox
from tkinter import filedialog, messagebox

import ttkbootstrap as ttk
from ttkbootstrap.constants import *

import json
import os

from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm



HONEY_DATA_FILE = "honey_data.json"
YEAST_DATA_FILE = "yeast_data.json"

def load_data(file_path, default_data):
    if os.path.exists(file_path):
        with open(file_path, "r") as f:
            return json.load(f)
    else:
        return default_data

def save_data(file_path, data):
    with open(file_path, "w") as f:
        json.dump(data, f, indent=4)

honey_data = load_data(HONEY_DATA_FILE, [
    {"Name": "Runny Honey Lidl", "Sugar": "79.7%", "Density": "1376.4 kg/m³", "Cost": "0.3"},
    {"Name": "Clear Honey Lidl", "Sugar": "79.7%", "Density": "1353.5 kg/m³", "Cost": "0.28"},
])

yeast_data = load_data(YEAST_DATA_FILE, [
    {"Name": "Lalvin 71B", "N Requirement": "Medium"},
    {"Name": "EC-1118", "N Requirement": "Low"},
])

def clear_screen(root):
    for widget in root.winfo_children():
        widget.destroy()

def show_home_screen(root):
    clear_screen(root)
    ttk.Label(root, text="Welcome to Mead Helper!", font=("Arial", 16)).pack(pady=20)
    ttk.Button(root, text="Mead Recipe", command=lambda: show_mead_recipe_screen(root)).pack(pady=20)
    ttk.Button(root, text="ABV measurement", command=lambda: show_ABV_calculation_screen(root)).pack(pady=20)
    ttk.Button(root, text="Backsweetening", command=lambda: show_backsweetening_calculation_screen(root)).pack(pady=20)
    ttk.Button(root, text="pH Adjustment", command=lambda: show_ph_adjustment_screen(root)).pack(pady=20)
    ttk.Button(root, text="Fermentation tracking", command=lambda: show_fermentation_tracking_screen(root)).pack(pady=20)
    root.iconbitmap("mead_calculation_icon.ico")

def show_mead_recipe_screen(root):
    clear_screen(root)

    top_bar = ttk.Frame(root)
    top_bar.pack(fill='x', padx=10, pady=5)

    ttk.Button(top_bar, text="Back", command=lambda: show_home_screen(root), bootstyle="success large").pack(side='left')
    ttk.Button(top_bar, text="Calculate", command=lambda: show_output_window(root), bootstyle="success large").pack(side='right')

    main_frame = ttk.Frame(root)
    main_frame.pack(fill='both', expand=True, padx=10, pady=10)

    # --- Recipe Form ---
    entry_group = ttk.LabelFrame(main_frame, text="Mead Recipe Form", padding=10)
    entry_group.grid(row=0, column=0, sticky='n', padx=10)

    root.batch_size_var = tk.StringVar()
    root.fg_var = tk.StringVar()
    root.abv_var = tk.StringVar()
    root.honey_type_var = tk.StringVar()
    root.yeast_type_var = tk.StringVar()
    root.use_fruit_var = tk.BooleanVar()
    root.fruit_type_var = tk.StringVar()

    ttk.Label(entry_group, text="Batch Size (L):").pack(anchor='w')
    ttk.Entry(entry_group, textvariable=root.batch_size_var).pack(fill='x', pady=4)

    ttk.Label(entry_group, text="Final Gravity (FG):").pack(anchor='w')
    ttk.Entry(entry_group, textvariable=root.fg_var).pack(fill='x', pady=4)

    ttk.Label(entry_group, text="Target ABV (%):").pack(anchor='w')
    ttk.Entry(entry_group, textvariable=root.abv_var).pack(fill='x', pady=4)

    ttk.Label(entry_group, text="Honey Type:").pack(anchor='w')
    honey_names = [h["Name"] for h in honey_data]
    honey_combo = ttk.Combobox(entry_group, textvariable=root.honey_type_var, values=honey_names, state="readonly")
    honey_combo.pack(fill='x', pady=4)

    ttk.Label(entry_group, text="Yeast Type:").pack(anchor='w')
    yeast_names = [y["Name"] for y in yeast_data]
    ttk.Combobox(entry_group, textvariable=root.yeast_type_var, values=yeast_names, state="readonly").pack(fill='x', pady=4)

    def toggle_fruit_entry():
        fruit_entry.config(state="normal" if root.use_fruit_var.get() else "disabled")

    ttk.Checkbutton(entry_group, text="Use Fruit?", variable=root.use_fruit_var, command=toggle_fruit_entry).pack(anchor='w', pady=(15, 0))
    fruit_entry = ttk.Entry(entry_group, textvariable=root.fruit_type_var, state="disabled")
    fruit_entry.pack(fill='x', pady=(0, 10))

    # --- Honey Databank with Headers ---
    honey_group = ttk.LabelFrame(main_frame, text="Honey Data", padding=10)
    honey_group.grid(row=0, column=1, sticky='n', padx=10)

    honey_canvas = tk.Canvas(honey_group, height=200)
    honey_scrollbar = ttk.Scrollbar(honey_group, orient="vertical", command=honey_canvas.yview)
    honey_frame = ttk.Frame(honey_canvas)

    honey_canvas.create_window((0, 0), window=honey_frame, anchor="nw")
    honey_canvas.configure(yscrollcommand=honey_scrollbar.set)
    honey_canvas.pack(side="left", fill="both", expand=True)
    honey_scrollbar.pack(side="right", fill="y")

    def refresh_honey():
        for widget in honey_frame.winfo_children():
            widget.destroy()
        headers = ["Name", "Sugar", "Density", "Cost", "Delete"]
        for c, h in enumerate(headers):
            ttk.Label(honey_frame, text=h, font=("Arial", 10, "bold")).grid(row=0, column=c, padx=5)
        for i, h in enumerate(honey_data):
            ttk.Label(honey_frame, text=h["Name"]).grid(row=i+1, column=0)
            ttk.Label(honey_frame, text=h["Sugar"]).grid(row=i+1, column=1)
            ttk.Label(honey_frame, text=h["Density"]).grid(row=i+1, column=2)
            ttk.Label(honey_frame, text=h.get("Cost", "N/A")).grid(row=i+1, column=3)
            ttk.Button(honey_frame, text="Delete", command=lambda i=i: delete_honey(i)).grid(row=i+1, column=4)

    def delete_honey(i):
        del honey_data[i]
        save_data(HONEY_DATA_FILE, honey_data)
        refresh_honey()

    refresh_honey()
    # --- Add Honey Form ---
    honey_form = ttk.LabelFrame(main_frame, text="Add New Honey", padding=10)
    honey_form.grid(row=1, column=1, sticky='n', padx=10)

    hname = tk.StringVar()
    hsugar = tk.StringVar()
    hdensity = tk.StringVar()
    hcost = tk.StringVar()

    ttk.Label(honey_form, text="Honey name:").pack(anchor='w')
    ttk.Entry(honey_form, textvariable=hname).pack(fill='x', pady=2)

    ttk.Label(honey_form, text="Sugar (%):").pack(anchor='w')
    ttk.Entry(honey_form, textvariable=hsugar).pack(fill='x', pady=2)

    ttk.Label(honey_form, text="Density (kg/m³):").pack(anchor='w')
    ttk.Entry(honey_form, textvariable=hdensity).pack(fill='x', pady=2)

    ttk.Label(honey_form, text="Cost (£/100g):").pack(anchor='w')
    ttk.Entry(honey_form, textvariable=hcost).pack(fill='x', pady=2)

    ttk.Button(honey_form, text="Add to Databank", command=lambda: add_honey(hname, hsugar, hdensity, hcost)).pack(pady=4)




    def add_honey(name, sugar, density, cost):
        try:
            honey_data.append({
                "Name": name.get(),
                "Sugar": f"{float(sugar.get())}%",
                "Density": f"{float(density.get())} kg/m³",
                "Cost": f"{float(cost.get()):.2f}"
            })
            save_data(HONEY_DATA_FILE, honey_data)
            refresh_honey()
            # Update dropdown values
            honey_combo["values"] = [h["Name"] for h in honey_data]
            name.set("")
            sugar.set("")
            density.set("")
            cost.set("")
        except Exception as e:
            tk.messagebox.showerror("Input Error", f"Please enter valid numbers.\n\n{e}")

    # --- Yeast Databank with Headers ---
    yeast_group = ttk.LabelFrame(main_frame, text="Yeast Data", padding=10)
    yeast_group.grid(row=0, column=2, sticky='n', padx=10)

    yeast_canvas = tk.Canvas(yeast_group, height=200)
    yeast_scrollbar = ttk.Scrollbar(yeast_group, orient="vertical", command=yeast_canvas.yview)
    yeast_frame = ttk.Frame(yeast_canvas)

    yeast_canvas.create_window((0, 0), window=yeast_frame, anchor="nw")
    yeast_canvas.configure(yscrollcommand=yeast_scrollbar.set)
    yeast_canvas.pack(side="left", fill="both", expand=True)
    yeast_scrollbar.pack(side="right", fill="y")

    def refresh_yeast():
        for widget in yeast_frame.winfo_children():
            widget.destroy()
        headers = ["Name", "N Requirement", "Delete"]
        for c, h in enumerate(headers):
            ttk.Label(yeast_frame, text=h, font=("Arial", 10, "bold")).grid(row=0, column=c, padx=5)
        for i, y in enumerate(yeast_data):
            ttk.Label(yeast_frame, text=y["Name"]).grid(row=i+1, column=0)
            ttk.Label(yeast_frame, text=y["N Requirement"]).grid(row=i+1, column=1)
            ttk.Button(yeast_frame, text="Delete", command=lambda i=i: delete_yeast(i)).grid(row=i+1, column=2)

    def delete_yeast(i):
        del yeast_data[i]
        save_data(YEAST_DATA_FILE, yeast_data)
        refresh_yeast()

    refresh_yeast()

    # --- Add Yeast Form ---
    yeast_form = ttk.LabelFrame(main_frame, text="Add New Yeast", padding=10)
    yeast_form.grid(row=1, column=2, sticky='n', padx=10)

    yname = tk.StringVar()
    ynreq = tk.StringVar()

    ttk.Label(yeast_form, text="Yeast name:").pack(anchor='w')
    ttk.Entry(yeast_form, textvariable=yname).pack(fill='x', pady=2)
    ttk.Label(yeast_form, text="N requirement:").pack(anchor='w')
    ttk.Combobox(yeast_form, textvariable=ynreq, values=["Low", "Medium", "High"], state="readonly").pack(fill='x', pady=2)
    ttk.Button(yeast_form, text="Add to Databank", command=lambda: add_yeast(yname, ynreq)).pack(pady=4)

    def add_yeast(name, nreq):
        if name.get().strip():
            yeast_data.append({"Name": name.get(), "N Requirement": nreq.get()})
            save_data(YEAST_DATA_FILE, yeast_data)
            refresh_yeast()
            name.set("")
            nreq.set("")
def show_output_window(root):
    try:
        volume_l = float(root.batch_size_var.get())
        final_gravity_desired = float(root.fg_var.get())
        ABV_desired = float(root.abv_var.get())
        honey_name = root.honey_type_var.get().strip()
        yeast_name = root.yeast_type_var.get().strip()
        fruit_used = root.use_fruit_var.get()
        fruit_type = root.fruit_type_var.get().strip()

        honey_match = next((h for h in honey_data if h["Name"].lower() == honey_name.lower()), None)
        sugar_conc = float(honey_match["Sugar"].replace('%', '').strip()) if honey_match else 0
        density = float(honey_match["Density"].replace('kg/m³', '').strip()) if honey_match else 1
        cost_per_100g = float(honey_match["Cost"]) if honey_match and "Cost" in honey_match else 0.3

        yeast_match = next((y for y in yeast_data if y["Name"].lower() == yeast_name.lower()), None)
        n_req = yeast_match["N Requirement"] if yeast_match else "Unknown"

        # Constants
        F_sp = 0.0128
        Y_xs = 0.1
        MW_CO2 = 44.01
        MW_eth = 46.069
        rho_eth = 789.45
        fraction_fermentable = 0.925

        # Gravity and sugar calculations
        starting_gravity = ((ABV_desired * 0.79) / (100 * 1.05) + 1)
        mass_ethanol = ((volume_l / 1000) * rho_eth * ABV_desired) / 100
        total_sugar_needed = ((1 / (1 - Y_xs)) * (mass_ethanol * (1 + (MW_CO2 / MW_eth)) +
                     F_sp * volume_l) * 1000)
        test_mass_honey = (((1 / (1 - Y_xs)) * (mass_ethanol * (1 + (MW_CO2 / MW_eth)) +
                              F_sp * volume_l) * 1000) / (sugar_conc / 100)) / fraction_fermentable if sugar_conc else 0

        total_honey_kg = test_mass_honey / 1000 if sugar_conc else 0
        volume_honey = total_honey_kg / (density / 1000) if density else 0
        cost = (test_mass_honey / 100) * cost_per_100g if sugar_conc else 0

        # Brix estimate
        brix = (
            182.46007 * starting_gravity**3
            - 775.68212 * starting_gravity**2
            + 1262.7794 * starting_gravity
            - 669.56218
        )

        # Fermaid-O calculation (ABV must be <= 14%)
        if ABV_desired <= 14:
            nitrogen_factors = {"Low": 0.75, "Medium": 0.9, "High": 1.25}
            N_req_val = nitrogen_factors.get(n_req, 0.9)
            volume_us_gallons = volume_l / 3.78541
            m_fermaid_o = ((brix * 10) * N_req_val * volume_us_gallons) / 50
        else:
            m_fermaid_o = None

        # Back-sweetening calculation
        imaginary_ABV_for_desired_final_sweetness = (1.05 / 0.79) * ((final_gravity_desired - 1) / 1) * 100
        mass_ethanol_sweetening = ((volume_l / 1000) * rho_eth * imaginary_ABV_for_desired_final_sweetness) / 100
        mass_pure_sugar_needed_for_sweetening = ((1 / (1 - Y_xs)) * (mass_ethanol_sweetening * (1 + (MW_CO2 / MW_eth)) +
                     F_sp * volume_l) * 1000)
        mass_honey_needed_for_sweetening = (((1 / (1 - Y_xs)) * (mass_ethanol_sweetening * (1 + (MW_CO2 / MW_eth)) +
                              F_sp * volume_l) * 1000) / (sugar_conc / 100)) / fraction_fermentable if sugar_conc else 0

        # Output
        output = tk.Toplevel(root)
        output.title("Calculation Output")
        output.geometry("1500x800")
        ttk.Label(output, text="Calculation Results", font=("Arial", 14)).pack(pady=10)

        result_text = f"""
Honey type: {honey_name or 'Unknown'}
Yeast: {yeast_name or 'Unknown'} (N Requirement: {n_req})
Fruit used: {'Yes - ' + fruit_type if fruit_used else 'No'}

Desired ABV: {ABV_desired:.1f}%
Final gravity target: {final_gravity_desired:.3f}
Starting gravity estimate: {starting_gravity:.3f}
Brix estimate: {brix:.1f}

Total pure sugar needed: {total_sugar_needed:.1f} g
Honey required: {total_honey_kg * 1000:.2f} g
Estimated cost: £{cost:.2f}
Volume of water to add: {volume_l - volume_honey:.2f} L
"""

        if m_fermaid_o is not None:
            result_text += f"""
--- Fermaid-O Nutrient Addition ---
Fermaid-O Required: {m_fermaid_o/4:.2f} g on day 0
Fermaid-O Required: {m_fermaid_o/4:.2f} g on day 1
Fermaid-O Required: {m_fermaid_o/4:.2f} g on day 2
Fermaid-O Required: {m_fermaid_o/4:.2f} g on day 3
This is all the nutrient that is required!
"""

        result_text += f"""
--- Honey for Back-Sweetening ---
Back-Sweetening Target FG: {final_gravity_desired:.3f}
Total pure sugar needed: {mass_pure_sugar_needed_for_sweetening:.2f} g
Honey Mass Required: {mass_honey_needed_for_sweetening:.2f} g
"""

        text_widget = tk.Text(output, wrap="word")
        text_widget.insert("1.0", result_text.strip())
        text_widget.config(state="disabled")
        text_widget.pack(fill="both", expand=True, padx=10)

        def save_as_pdf():
            file_path = filedialog.asksaveasfilename(
                defaultextension = ".pdf",
                filetypes = [("PDF files", "*.pdf")],
                title = "Save as PDF"
            )
            if file_path:
                try:
                    from reportlab.pdfgen import canvas
                    from reportlab.lib.pagesizes import A4
                    from reportlab.lib.units import mm

                    c = canvas.Canvas(file_path, pagesize=A4)
                    width, height = A4
                    lines = result_text.strip().split('\n')
                    y = height - 20 * mm

                    for line in lines:
                        c.drawString(20 * mm, y, line)
                        y -= 10 * mm
                        if y < 20 * mm:
                            c.showPage()
                            y = height - 20 * mm
                    
                    c.save()
                    messagebox.showinfo("Success", f"PDF saved successfully:\n{file_path}")
                except Exception as e:
                    messagebox.showerror("Error", f"Failed to save PDF:\n{e}")

        
        ttk.Button(output, text="Save as PDF", command=save_as_pdf).pack(pady=10)
        ttk.Button(output, text="Close", command=output.destroy).pack(pady=10)

    
    except Exception as e:
        tk.messagebox.showerror("Error", f"Calculation failed:\n{e}")







def show_ABV_calculation_screen(root):
    clear_screen(root)

    ttk.Label(root, text="ABV calculation", font=("Arial", 16)).pack(pady=10)

    container = ttk.Frame(root)
    container.pack(pady=10, padx=20, anchor="center")

    starting_gravity = tk.StringVar()
    final_gravity = tk.StringVar()

    form = ttk.Frame(container)
    form.grid(row=0, column=0, sticky="n", padx=(0, 20))

    ttk.Label(form, text="Starting Gravity:").pack(anchor='w', pady=(0, 2))
    ttk.Entry(form, textvariable=starting_gravity, width=30).pack(fill='x', pady=(0, 10))

    ttk.Label(form, text="Final Gravity:").pack(anchor='w', pady=(0, 2))
    ttk.Entry(form, textvariable=final_gravity, width=30).pack(fill='x', pady=(0, 10))

    btn_frame = ttk.Frame(root)
    btn_frame.pack(pady=10)

    def calculate_ABV():
        try:
            starting_gravity_calc = float(starting_gravity.get())
            final_gravity_calc = float(final_gravity.get())

            ABV = (1.05/0.79) * ((starting_gravity_calc - final_gravity_calc)/(final_gravity_calc)) * 100

            result = f"""
Starting Gravity: {starting_gravity_calc}
Final Gravity: {final_gravity_calc}

%ABV calculated: {ABV:.1f} %
"""

            output = tk.Toplevel(root)
            output.title("ABV Calculation Result")
            output.geometry("1400x1300")
            ttk.Label(output, text="ABV Calculation", font=("Arial", 14)).pack(pady=10)
            text = tk.Text(output, wrap="word")
            text.insert("1.0", result.strip())
            text.config(state="disabled")
            text.pack(fill="both", expand=True, padx=10)
            ttk.Button(output, text="Close", command=output.destroy).pack(pady=10)

        except Exception as e:
            tk.messagebox.showerror("Error", f"Invalid input:\n{e}")

    ttk.Button(btn_frame, text="Calculate", width=20, command=calculate_ABV, bootstyle="success large").pack(pady=5)
    ttk.Button(btn_frame, text="Back", width=20, command=lambda: show_home_screen(root), bootstyle="success large").pack()

def show_ph_adjustment_screen(root):
    clear_screen(root)

    ttk.Label(root, text="pH Adjustment", font=("Arial", 16)).pack(pady=10)

    # Use center container
    container = ttk.Frame(root)
    container.pack(pady=10, padx=20, anchor="center")

    current_ph = tk.StringVar()
    target_ph = tk.StringVar()
    volume_l = tk.StringVar()

    # Form column
    form = ttk.Frame(container)
    form.grid(row=0, column=0, sticky="n", padx=(0, 20))

    ttk.Label(form, text="Current pH:").pack(anchor='w', pady=(0, 2))
    ttk.Entry(form, textvariable=current_ph, width=30).pack(fill='x', pady=(0, 10))

    ttk.Label(form, text="Target pH:").pack(anchor='w', pady=(0, 2))
    ttk.Entry(form, textvariable=target_ph, width=30).pack(fill='x', pady=(0, 10))

    ttk.Label(form, text="Volume (L):").pack(anchor='w', pady=(0, 2))
    ttk.Entry(form, textvariable=volume_l, width=30).pack(fill='x', pady=(0, 10))

    # Notes side panel
    notes_frame = ttk.LabelFrame(container, text="Info / Notes")
    notes_frame.grid(row=0, column=1, sticky="n")

    notes_box = tk.Text(notes_frame, width=40, height=10, wrap="word")
    notes_box.pack(fill='both', expand=True)

    notes_box.insert("1.0", "This calculation is based on calcium carbonate, any other pH adjusting compound will use a different calculation.\n")

    # Bottom buttons
    btn_frame = ttk.Frame(root)
    btn_frame.pack(pady=10)

    def calculate_ph_adjustment():
        try:
            pH_initial = float(current_ph.get())
            pH_target = float(target_ph.get())
            V = float(volume_l.get())

            H_initial = 10 ** (-pH_initial)
            H_target = 10 ** (-pH_target)
            delta_H = (H_initial - H_target) * V * 1000
            mol_CaCO3 = delta_H / 2
            mass_CaCO3 = mol_CaCO3 * 100.09

            result = f"""
Initial pH: {pH_initial}
Target pH: {pH_target}
Volume: {V} L

[H⁺] Initial: {H_initial:.2e} mol/L
[H⁺] Target: {H_target:.2e} mol/L
Δ[H⁺] Total: {delta_H:.2e} mol

Moles CaCO₃ needed: {mol_CaCO3:.4f} mol
Mass CaCO₃ required: {mass_CaCO3:.4f} g
"""

            output = tk.Toplevel(root)
            output.title("pH Adjustment Result")
            output.geometry("1400x1300")
            ttk.Label(output, text="CaCO₃ Requirement", font=("Arial", 14)).pack(pady=10)
            text = tk.Text(output, wrap="word")
            text.insert("1.0", result.strip())
            text.config(state="disabled")
            text.pack(fill="both", expand=True, padx=10)
            ttk.Button(output, text="Close", command=output.destroy).pack(pady=10)

        except Exception as e:
            tk.messagebox.showerror("Error", f"Invalid input:\n{e}")

    ttk.Button(btn_frame, text="Calculate", width=20, command=calculate_ph_adjustment, bootstyle="success large").pack(pady=5)
    ttk.Button(btn_frame, text="Back", width=20, command=lambda: show_home_screen(root), bootstyle="success large").pack()


def show_fermentation_tracking_screen(root):

    clear_screen(root)

    ttk.Label(root, text="Fermentation Tracking", font=("Arial", 16)).pack(pady=10)

    container = ttk.Frame(root)
    container.pack(pady=10, padx=20, fill="both", expand=True)

    # Left Side - Entry Table
    left_frame = ttk.Frame(container)
    left_frame.pack(side="left", padx=20, pady=10)

    ttk.Label(left_frame, text="Time (days)").grid(row=0, column=0, padx=5, pady=5)
    ttk.Label(left_frame, text="SG").grid(row=0, column=1, padx=5, pady=5)

    time_entries = []
    sg_entries = []

    for i in range(6):
        time_entry = ttk.Entry(left_frame, width=10)
        sg_entry = ttk.Entry(left_frame, width=10)
        time_entry.grid(row=i+1, column=0, padx=5, pady=2)
        sg_entry.grid(row=i+1, column=1, padx=5, pady=2)
        time_entries.append(time_entry)
        sg_entries.append(sg_entry)

    # Right Side - Graph Area
    right_frame = ttk.Frame(container)
    right_frame.pack(side="right", padx=20, pady=30, fill="both", expand=True)

    fig = Figure(figsize=(6, 4), dpi=100)
    ax1 = fig.add_subplot(211)  # Top graph
    ax2 = fig.add_subplot(212)  # Bottom graph

    fig.tight_layout(pad=2.0)

    canvas = FigureCanvasTkAgg(fig, master=right_frame)
    canvas.get_tk_widget().pack(fill="both", expand=True)


    def calculate_graphs():
        from scipy.optimize import differential_evolution

        try:
            times = []
            sgs = []
            for i in range(6):
                try:
                    t = float(time_entries[i].get())
                    sg = float(sg_entries[i].get())
                    times.append(t)
                    sgs.append(sg)
                except ValueError:
                    continue

            if len(times) < 2:
                tk.messagebox.showerror("Input Error", "Please enter at least 2 valid data points.")
                return

            # Sort inputs
            paired = sorted(zip(times, sgs))
            times, sgs = zip(*paired)

            SG_min = 1.000
            SG_max = max(sgs)
            SG0 = float(sgs[0])
            V_mead = float(volume_of_mead.get())
            X0 = float(mass_of_yeast.get()) / V_mead
            YXS = 0.1  # <-- FIXED yield coefficient

            # Constants
            F_SP = 0.0128
            rho_eth = 789.45
            MW_CO2 = 44.01
            MW_eth = 46.069

            # Initial sugar concentration (fixed)
            S0 = V_mead * (((SG0 - 1) * ((1.05 / 0.79) * rho_eth) * (1 + (MW_CO2 / MW_eth)) / (1 - YXS)) + F_SP)


            t_max = max(times + (float(predict_time_entry.get()),))
            t_fit = np.linspace(0, t_max + 2, 600)

            # --- Logistic Fit ---
            def logistic_eqns(x):
                k, t0 = x
                return [
                    SG_min + (SG_max - SG_min) / (1 + np.exp(-k * (times[0] - t0))) - sgs[0],
                    SG_min + (SG_max - SG_min) / (1 + np.exp(-k * (times[-1] - t0))) - sgs[-1]
                ]

            try:
                k, t0 = fsolve(logistic_eqns, [-2.0, times[-1] / 2])
                sg_fit = SG_min + (SG_max - SG_min) / (1 + np.exp(-k * (t_fit - t0)))

                ax1.clear()
                ax1.plot(t_fit, sg_fit, 'b-', label='Logistic Fit')
                ax1.scatter(times, sgs, color='red', label='Measured')
                ax1.set_title("Graph 1: Logistic SG Fit")
                ax1.set_xlabel("Time (days)")
                ax1.set_ylabel("SG")
                ax1.set_ylim(0.98, SG0*1.05)
                ax1.grid(True)
                ax1.legend()
            except Exception as e:
                ax1.clear()
                ax1.text(0.1, 0.5, f"Logistic fit failed:\n{e}", color='red')

            # --- Monod Fit (only mu_max and Ks) ---
            def simulate_monod(mu_max, Ks, t_eval):
                def dXdt(t, y):
                    X, S = y
                    S = max(S, 1e-6)
                    mu = mu_max * S / (Ks + S + 1e-8)
                    dX = mu * X
                    dS = -dX / YXS if S > 1e-6 else 0
                    return [dX, dS]

                try:
                    sol = solve_ivp(dXdt, [0, max(t_eval)], [X0, S0], t_eval=t_eval, max_step=0.1)
                    S = np.maximum(sol.y[1], 0)
                    SG = 1 + (1 - YXS) * (S / V_mead - F_SP) / (((1.05 / 0.79) * rho_eth) * (1 + (MW_CO2 / MW_eth)))
                    SG = np.maximum(SG, SG_min)
                    return SG, sol.t
                except:
                    return np.full_like(t_eval, np.nan), t_eval

            def objective(params):
                mu_max, Ks = params
                if not (0.001 < mu_max < 5 and 0.01 < Ks < 50):
                    return np.inf
                SG_pred, _ = simulate_monod(mu_max, Ks, np.array(times))
                if np.any(np.isnan(SG_pred)):
                    return np.inf
                return np.sum((SG_pred - np.array(sgs)) ** 2)

            result = differential_evolution(objective, bounds=[(0.001, 5), (0.01, 50)],
                                        strategy='best1bin', maxiter=500, polish=True)

            if result.success:
                mu_opt, Ks_opt = result.x
                SG_fit, t_vals = simulate_monod(mu_opt, Ks_opt, t_fit)

                ax2.clear()
                ax2.plot(t_vals, SG_fit, color='green', label='Monod Fit')
                ax2.scatter(times, sgs, color='red', label='Measured')
                ax2.set_title("Graph 2: Fitted Monod SG Model")
                ax2.set_xlabel("Time (days)")
                ax2.set_ylabel("SG")
                ax2.set_ylim(0.98, SG0*1.05)
                ax2.grid(True)
                ax2.legend()

                # Predict SG at a specified time
                try:
                    predict_time = float(predict_time_entry.get())
                    SG_pred, _ = simulate_monod(mu_opt, Ks_opt, [predict_time])
                    if np.isnan(SG_pred[-1]):
                        raise ValueError
                    sg_monod_label.config(text=f"Monod SG: {SG_pred[-1]:.4f}")
                except:
                    sg_monod_label.config(text="Monod SG: —")
            else:
                ax2.clear()
                ax2.text(0.1, 0.5, "Monod fit failed", color='red')
                sg_monod_label.config(text="Monod SG: —")

            # Predict SG from Logistic model
            try:
                predict_time = float(predict_time_entry.get())
                sg_log = SG_min + (SG_max - SG_min) / (1 + np.exp(-k * (predict_time - t0)))
                sg_logistic_label.config(text=f"Logistic SG: {sg_log:.4f}")
            except:
                sg_logistic_label.config(text="Logistic SG: —")

            canvas.draw()

        except Exception as e:
            tk.messagebox.showerror("Error", str(e))

    # Bottom Buttons
    ttk.Button(left_frame, text="Calculate", command=calculate_graphs, bootstyle="success large").grid(row=7, column=0, columnspan=2, pady=10)

    # inital concentration of yeast

    ttk.Label(left_frame, text="Mass of yeast added (g)").grid(row=8, column=0, pady=5)
    mass_of_yeast = ttk.Entry(left_frame, width=10)
    mass_of_yeast.grid(row=8, column=1, padx=5)
    ttk.Label(left_frame, text="Volume of mead being made (L)").grid(row=10, column=0)
    volume_of_mead = ttk.Entry(left_frame, width=10)
    volume_of_mead.grid(row=10, column=1, pady=10)

        # Time prediction input
    ttk.Label(left_frame, text="Predict SG at time:").grid(row=11, column=0, padx=5, pady=(15, 2), sticky="w")
    predict_time_entry = ttk.Entry(left_frame, width=10)
    predict_time_entry.grid(row=11, column=1, padx=5, pady=(15, 2))

    sg_logistic_label = ttk.Label(left_frame, text="Logistic SG: —")
    sg_logistic_label.grid(row=12, column=0, columnspan=2, sticky="w", padx=5)

    sg_monod_label = ttk.Label(left_frame, text="Monod SG: —")
    sg_monod_label.grid(row=13, column=0, columnspan=2, sticky="w", padx=5)


    btn_frame = ttk.Frame(root)
    btn_frame.pack(pady=10)
    ttk.Button(btn_frame, text="Back", width=20, command=lambda: show_home_screen(root), bootstyle="success large").pack()




def show_backsweetening_calculation_screen(root):
    clear_screen(root)

    ttk.Label(root, text="Back-sweetening Calculation", font=("Arial", 16)).pack(pady=10)

    # Use center container
    container = ttk.Frame(root)
    container.pack(pady=10, padx=20, anchor="center")

    final_gravity = tk.StringVar()
    desired_gravity = tk.StringVar()
    volume_l = tk.StringVar()

    # Form column
    form = ttk.Frame(container)
    form.grid(row=0, column=0, sticky="n", padx=(0, 20))

    ttk.Label(form, text="Final Gravity reading:").pack(anchor='w', pady=(0, 2))
    ttk.Entry(form, textvariable=final_gravity, width=30).pack(fill='x', pady=(0, 10))

    ttk.Label(form, text="Target Gravity:").pack(anchor='w', pady=(0, 2))
    ttk.Entry(form, textvariable=desired_gravity, width=30).pack(fill='x', pady=(0, 10))

    ttk.Label(form, text="Volume of mead (L):").pack(anchor='w', pady=(0, 2))
    ttk.Entry(form, textvariable=volume_l, width=30).pack(fill='x', pady=(0, 10))

    # Bottom buttons
    btn_frame = ttk.Frame(root)
    btn_frame.pack(pady=10)

    def calculate_backsweetening():
        try:
            final_gravity_reading = float(final_gravity.get())
            target_gravity = float(desired_gravity.get())
            V = float(volume_l.get())

            mass_sugar_needed = 
            mass_honey_needed = 

            result = f"""
Final Gravity: {final_gravity_reading}
Fina Gravity desired: {target_gravity}
Volume: {V} L

Mass of sugar needed: {mass_sugar_needed:.2e} g
Mass of honey needed: {mass_honey_needed:.2e} g
"""

            output = tk.Toplevel(root)
            output.title("Back-sweetening Result")
            output.geometry("1400x1300")
            ttk.Label(output, text="Sugar/Honey requirement", font=("Arial", 14)).pack(pady=10)
            text = tk.Text(output, wrap="word")
            text.insert("1.0", result.strip())
            text.config(state="disabled")
            text.pack(fill="both", expand=True, padx=10)
            ttk.Button(output, text="Close", command=output.destroy).pack(pady=10)

        except Exception as e:
            tk.messagebox.showerror("Error", f"Invalid input:\n{e}")

    ttk.Button(btn_frame, text="Calculate", width=20, command=calculate_backsweetening, bootstyle="success large").pack(pady=5)
    ttk.Button(btn_frame, text="Back", width=20, command=lambda: show_home_screen(root), bootstyle="success large").pack()

def main():
    root = tk.Tk()
    root.title("Mead Assistant")

    # Start full-screen
    root.attributes("-fullscreen", True)

    # Function to exit fullscreen and resize + center
    def exit_fullscreen(event=None):
        root.attributes("-fullscreen", False)
        w, h = 1300, 730
        screen_width = root.winfo_screenwidth()
        screen_height = root.winfo_screenheight()
        x = (screen_width // 2) - (w // 2)
        y = (screen_height // 2) - (h // 2)
        root.geometry(f"{w}x{h}+{x}+{y}")

    # Bind Esc key to exit fullscreen
    root.bind("<Escape>", exit_fullscreen)

    show_home_screen(root)
    root.mainloop()


if __name__ == "__main__":
    main()

# python -m PyInstaller --noconfirm --onefile --windowed --icon="mead_calculation_icon.ico" --collect-submodules=numpy --collect-submodules=scipy "Homebrewing calculator.py"