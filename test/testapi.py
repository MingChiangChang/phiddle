import requests
import numpy as np
import matplotlib.pyplot as plt


data = {}
data["std_noise"] = 0.05 
data["mean_θ"] = [1.0, 0.5, 0.2]
data["std_θ"] = [0.05, 0.05, 0.05]
data["maxiter"] = 128
data["regularization"] = True 
data["method"] = "LM"
data["objective"] = "LS"
data["optimization_mode"] = "Simple"
data["em_loop_num"] = 5
data["λ"] = 1. 
data["verbose"] = False
data["tol"] = 1E-6
data["data"] = np.random.rand(1000).tolist()
data["q"] = np.linspace(10, 80, 1000).tolist()
data["depth"] = 1
data["k"] = 1
data["amorphous"] = False
data["background"] = False
data["background_length"] = 5.
data["background_option"] = "MCBL"
data["csv_file"] = "/Users/ming/Desktop/Data/cifsssss/SnOx/test.csv"
data["names"] = ["SnO2_P42/mnm", "Sn3O4_P1"]
data["phase_name"] = "SnO2_P42/mnm"

response = requests.put("http://127.0.0.1:8080/set_csv", json=data)
response = requests.get("http://127.0.0.1:8080/phase_names")
print(response.json())
response = requests.get("http://127.0.0.1:8080/peak_info", json=data)
print(response.json())
response = requests.post("http://127.0.0.1:8080/label", json=data)
print(response.json())
# response = requests.post("http://127.0.0.1:8080/label_with_phase", json=data)
# response = requests.post("http://127.0.0.1:8080/fit_with_phase", json=data)


response = requests.get("http://127.0.0.1:8080/evaluate", json=data)
d = response.json()
plt.plot(data["q"], d)
plt.show()

response = requests.get("http://127.0.0.1:8080/phase_names")
phase_names = response.json()

for i in range(5):
    response = requests.get(f"http://127.0.0.1:8080/evaluate/{i}", json=data)
    d = response.json()
    
    plt.plot(data["q"], d)
    plt.title(phase_names[i])
    plt.show()


for i in range(5):
    data["phase_name"] = phase_names[i]
    response = requests.get(f"http://127.0.0.1:8080/evaluate", json=data)
    d = response.json()
    
    plt.plot(data["q"], d)
    plt.title(phase_names[i])
    plt.show()


for i in range(5):
    data["idx"] = i
    response = requests.get(f"http://127.0.0.1:8080/evaluate_fitted/{i}", json=data)
    print(response)
    d = response.json()
    response = requests.get(f"http://127.0.0.1:8080/get_fitted_phase_name/{i}", json=data)
    print(response)
    names = response.json()
    
    plt.plot(data["q"], d)
    plt.plot(data["q"], data["data"])
    plt.title(' '.join(names))
    plt.show()



