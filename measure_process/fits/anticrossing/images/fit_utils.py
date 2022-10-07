from dataclasses import dataclass
from typing import Dict, List

@dataclass
class ModelParameter:
    initial: float
    minimum: float
    maximum: float

class ModelParameters:
    def __init__(self, model_params:Dict[str,ModelParameter], fit_params:List[str]):
        self.model_params = model_params
        self.fit_params = fit_params

    def get_initial_values(self):
        return [self.model_params[p].initial for p in self.fit_params]

    def get_boundaries(self):
        return [(self.model_params[p].minimum, self.model_params[p].maximum)
                for p in self.fit_params]

    def get_model_values(self, fit_values:List[float]):
        values = {k:v.initial for k,v in self.model_params.items()}
        values.update(self.get_fitted_values(fit_values))
        return values

    def get_fitted_values(self, fit_values:List[float]):
        return {p:fit_values[i] for i,p in enumerate(self.fit_params)}
