# NewtonGrid

**Projeto:** Visualização Numérica de Campos Gravitacionais (VTK-ready)

## Descrição

O NewtonGrid é um núcleo numérico em C++ para calcular campos gravitacionais newtonianos em domínios 2D ou 3D regulares. O projeto permite simular múltiplas massas, validar propriedades físicas e exportar resultados para análise rápida ou visualização científica.

## Estrutura do Projeto

```
NewtonGrid/
├─ CMakeLists.txt            # Build configuration
├─ README.md                 # Este arquivo
├─ docs/                     # Documentação do projeto
├─ src/                      # Código fonte
│  ├─ main.cpp               # CLI entrypoint para experimentos
│  ├─ app/                   # Orquestra execução de experimentos
│  ├─ physics/               # Módulos de física (Mass, GravitationalField)
│  ├─ data/                  # Grid e FieldData
│  ├─ io/                    # Exportação CSV e VTI
│  ├─ utils/                 # Config, logger e utilitários matemáticos
│  └─ tests/                 # Testes unitários
├─ examples/                 # Configurações e outputs de exemplo
├─ extern/                   # Bibliotecas de terceiros (opcional)
└─ ci/                       # Integração contínua
```

## Objetivo

* Simular campos gravitacionais em 2D ou 3D.
* Fornecer resultados exportáveis para análise e visualização (CSV, VTI).
* Garantir modularidade e extensibilidade para diferentes experimentos e validações.

## Funcionamento

1. O usuário fornece uma configuração JSON definindo domínio, massas e parâmetros gravitacionais.
2. O `experiment_runner` lê a configuração e constrói o grid.
3. O campo gravitacional é calculado ponto a ponto.
4. Resultados são exportados para CSV e VTI.
5. Testes podem ser executados para validação numérica do cálculo.

