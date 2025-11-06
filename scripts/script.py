#!/usr/bin/env python3
"""
Network Propagation Script - HAB Tarea 2
=========================================

Este script implementa algoritmos de propagación en redes para priorización de genes:
1. GUILD (Random Walk with Restart)
2. DIAMOnD (DIseAse MOdule Detection)

Autor: Aissa Omar El Hammouti Chachoui
Uso: python tu_script.py --algorithm [guild|diamond|both] --network <file> --seeds <file> --output <file>
"""

import argparse
import sys
import pandas as pd
import networkx as nx
import numpy as np
from typing import List, Dict, Set, Tuple
import logging

# Configuración del logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class NetworkPropagation:
    """
    Clase para realizar propagación en redes usando diferentes algoritmos.
    
    Attributes:
        network_file (str): Ruta al archivo de red
        seeds_file (str): Ruta al archivo de genes semilla
        G (nx.Graph): Grafo de la red de interacciones
        seeds (set): Conjunto de genes semilla
    """
    
    def __init__(self, network_file: str, seeds_file: str):
        """
        Inicializa la clase con los archivos de red y semillas.
        
        Args:
            network_file: Archivo con la red de interacciones
            seeds_file: Archivo con los genes semilla
        """
        self.network_file = network_file
        self.seeds_file = seeds_file
        self.G = None
        self.seeds = set()
        
    def load_network(self, network_type: str = 'auto') -> nx.Graph:
        """
        Carga la red de interacciones desde el archivo.
        
        Soporta diferentes formatos:
        - GUILD format: node1 node2 weight
        - DIAMOnD format: node1 node2 weight
        - STRING format: node1 node2 combined_score
        
        Args:
            network_type: Tipo de formato ('guild', 'diamond', 'string', 'auto')
            
        Returns:
            nx.Graph: Grafo con las interacciones cargadas
        """
        logger.info(f"Cargando red desde {self.network_file}...")
        
        try:
            # Intentar cargar el archivo
            if network_type == 'auto':
                # Detectar formato automáticamente
                with open(self.network_file, 'r') as f:
                    first_line = f.readline().strip()
                    if '\t' in first_line:
                        sep = '\t'
                    else:
                        sep = ' '
            else:
                sep = '\t' if 'string' in network_type else ' '
            
            # Cargar datos
            df = pd.read_csv(
                self.network_file,
                sep=sep,
                header=None,
                names=['node1', 'node2', 'weight'],
                comment='#'
            )
            
            # Crear grafo no dirigido con pesos
            self.G = nx.Graph()
            
            for _, row in df.iterrows():
                node1 = str(row['node1']).strip()
                node2 = str(row['node2']).strip()
                weight = float(row['weight'])
                
                self.G.add_edge(node1, node2, weight=weight)
            
            logger.info(f"Red cargada: {self.G.number_of_nodes()} nodos, "
                       f"{self.G.number_of_edges()} aristas")
            
            return self.G
            
        except Exception as e:
            logger.error(f"Error al cargar la red: {e}")
            raise
    
    def load_seeds(self) -> Set[str]:
        """Carga genes semilla y los adapta automáticamente a la red."""
        logger.info(f"Cargando genes semilla desde {self.seeds_file}...")
    
        # Mapeo bidireccional nombre ↔ ID
        gene_map = {
            'ENO1': '2023', '2023': 'ENO1',
            'PGK1': '5230', '5230': 'PGK1',
            'HK2': '3099', '3099': 'HK2'
        }
        
        try:
            with open(self.seeds_file, 'r') as f:
                raw_seeds = [line.strip() for line in f if line.strip()]
            
            network_nodes = set(self.G.nodes())
            self.seeds = set()
            
            # Intentar cargar semillas con conversión automática
            for seed in raw_seeds:
                if seed in network_nodes:
                    self.seeds.add(seed)
                elif seed in gene_map and gene_map[seed] in network_nodes:
                    converted = gene_map[seed]
                    self.seeds.add(converted)
                    logger.info(f"Convertido: {seed} → {converted}")
            
            # Si no hay semillas, usar los 3 nodos con mayor grado
            if not self.seeds:
                logger.warning("Semillas no encontradas. Usando nodos de mayor grado...")
                degrees = dict(self.G.degree())
                top_nodes = sorted(degrees.items(), key=lambda x: x[1], reverse=True)[:3]
                self.seeds = {node for node, _ in top_nodes}
                logger.info(f"Usando como semillas: {self.seeds}")
            
            logger.info(f"Genes semilla cargados: {len(self.seeds)}")
            return self.seeds
            
        except Exception as e:
            logger.error(f"Error: {e}")
            raise
    
    def guild_random_walk_restart(self, alpha: float = 0.85, 
                                  max_iter: int = 100, 
                                  tol: float = 1e-6) -> Dict[str, float]:
        """
        Implementa el algoritmo GUILD usando Random Walk with Restart (RWR).
        
        El método RWR es uno de los más populares en GUILD. Simula un caminante
        aleatorio que en cada paso puede:
        1. Moverse a un nodo vecino (con probabilidad alpha)
        2. Regresar a un nodo semilla (con probabilidad 1-alpha)
        
        Matemáticamente: p(t+1) = alpha * W * p(t) + (1-alpha) * p(0)
        Donde:
        - p(t): vector de probabilidades en el tiempo t
        - W: matriz de transición normalizada
        - p(0): vector inicial (semillas tienen valor uniforme)
        
        Args:
            alpha: Parámetro de restart (típicamente 0.7-0.9)
            max_iter: Número máximo de iteraciones
            tol: Tolerancia para convergencia
            
        Returns:
            dict: Diccionario con scores para cada gen {gen: score}
        """
        logger.info(f"Ejecutando GUILD Random Walk with Restart (alpha={alpha})...")
        
        if not self.seeds:
            raise ValueError("No hay genes semilla cargados")
        
        # Obtener nodos y crear índices
        nodes = list(self.G.nodes())
        n_nodes = len(nodes)
        node_to_idx = {node: idx for idx, node in enumerate(nodes)}
        
        # Crear matriz de adyacencia ponderada
        A = nx.to_numpy_array(self.G, nodelist=nodes, weight='weight')
        
        # Normalizar por grados (matriz de transición)
        # W[i,j] = A[i,j] / sum_k(A[i,k])
        row_sums = A.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1  # Evitar división por cero
        W = A / row_sums
        
        # Vector inicial: distribución uniforme en semillas
        p0 = np.zeros(n_nodes)
        seed_indices = [node_to_idx[seed] for seed in self.seeds if seed in node_to_idx]
        
        if seed_indices:
            p0[seed_indices] = 1.0 / len(seed_indices)
        else:
            raise ValueError("Ninguna semilla encontrada en la red")
        
        # Iteración del Random Walk
        p = p0.copy()
        
        for iteration in range(max_iter):
            p_new = alpha * W.T @ p + (1 - alpha) * p0
            
            # Verificar convergencia
            diff = np.linalg.norm(p_new - p)
            
            if diff < tol:
                logger.info(f"Convergencia alcanzada en iteración {iteration + 1}")
                break
                
            p = p_new
        else:
            logger.warning(f"No se alcanzó convergencia después de {max_iter} iteraciones")
        
        # Crear diccionario de scores
        scores = {nodes[i]: float(p[i]) for i in range(n_nodes)}
        
        logger.info("GUILD completado")
        return scores
    
    def diamond_algorithm(self, max_genes: int = 200, 
                          alpha: float = 1.0) -> List[Tuple[str, int, float]]:
        """
        Implementa el algoritmo DIAMOnD (DIseAse MOdule Detection).
        
        DIAMOnD es un algoritmo iterativo que expande el conjunto de genes semilla:
        1. Comienza con las semillas
        2. En cada iteración, calcula para cada gen candidato su conectividad
           al módulo actual
        3. Añade el gen con mayor conectividad
        4. Repite hasta alcanzar el tamaño deseado
        
        La conectividad se mide como:
        k_s = suma de pesos de aristas entre el gen y el módulo
        
        Referencia: Ghiassian et al. (2015) PLOS Comput Biol
        
        Args:
            max_genes: Número máximo de genes a añadir
            alpha: Parámetro de penalización por grado (no usado en versión básica)
            
        Returns:
            list: Lista de tuplas (gen, iteración, score)
        """
        logger.info(f"Ejecutando DIAMOnD (max_genes={max_genes})...")
        
        if not self.seeds:
            raise ValueError("No hay genes semilla cargados")
        
        # Inicializar módulo con semillas
        module = set(self.seeds)
        all_nodes = set(self.G.nodes())
        candidates = all_nodes - module
        
        results = []
        
        # Añadir semillas a resultados con iteración 0
        for seed in self.seeds:
            results.append((seed, 0, float('inf')))
        
        # Iteraciones de DIAMOnD
        for iteration in range(1, max_genes + 1):
            if not candidates:
                logger.info(f"No quedan candidatos después de {iteration - 1} iteraciones")
                break
            
            # Calcular conectividad de cada candidato al módulo
            max_connectivity = -1
            best_candidate = None
            
            for candidate in candidates:
                # Sumar pesos de aristas entre el candidato y el módulo
                connectivity = 0
                for module_gene in module:
                    if self.G.has_edge(candidate, module_gene):
                        connectivity += self.G[candidate][module_gene].get('weight', 1.0)
                
                if connectivity > max_connectivity:
                    max_connectivity = connectivity
                    best_candidate = candidate
            
            # Añadir mejor candidato al módulo
            if best_candidate:
                module.add(best_candidate)
                candidates.remove(best_candidate)
                results.append((best_candidate, iteration, max_connectivity))
                
                if iteration % 50 == 0:
                    logger.info(f"Iteración {iteration}: {best_candidate} "
                              f"(conectividad={max_connectivity:.3f})")
            else:
                break
        
        logger.info(f"DIAMOnD completado: {len(results)} genes en el módulo")
        return results
    
    def save_results(self, results: Dict[str, any], output_file: str, 
                     algorithm: str) -> None:
        """
        Guarda los resultados en un archivo.
        
        Args:
            results: Resultados del algoritmo (dict o list)
            output_file: Ruta del archivo de salida
            algorithm: Nombre del algoritmo usado
        """
        logger.info(f"Guardando resultados en {output_file}...")
        
        try:
            if algorithm == 'guild':
                # Ordenar por score descendente
                sorted_results = sorted(results.items(), 
                                      key=lambda x: x[1], 
                                      reverse=True)
                
                df = pd.DataFrame(sorted_results, columns=['Gene', 'GUILD_Score'])
                df['Rank'] = range(1, len(df) + 1)
                df = df[['Rank', 'Gene', 'GUILD_Score']]
                
            elif algorithm == 'diamond':
                df = pd.DataFrame(results, 
                                columns=['Gene', 'DIAMOnD_Iteration', 'Connectivity'])
                df['Rank'] = range(1, len(df) + 1)
                df = df[['Rank', 'Gene', 'DIAMOnD_Iteration', 'Connectivity']]
            
            # Guardar a archivo
            df.to_csv(output_file, sep='\t', index=False)
            logger.info(f"Resultados guardados exitosamente ({len(df)} genes)")
            
            # Mostrar top 10
            logger.info(f"\nTop 10 genes:")
            print(df.head(10).to_string(index=False))
            
        except Exception as e:
            logger.error(f"Error al guardar resultados: {e}")
            raise


def main():
    """
    Función principal que maneja la interfaz de línea de comandos.
    """
    parser = argparse.ArgumentParser(
        description='Network Propagation usando GUILD y/o DIAMOnD',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos de uso:
  # Ejecutar solo GUILD
  python tu_script.py --algorithm guild --network data/network_guild.txt \\
      --seeds data/genes_seed.txt --output results/guild_results.txt
  
  # Ejecutar solo DIAMOnD
  python tu_script.py --algorithm diamond --network data/network_diamond.txt \\
      --seeds data/genes_seed.txt --output results/diamond_results.txt
  
  # Ejecutar ambos algoritmos
  python tu_script.py --algorithm both --network data/network_guild.txt \\
      --seeds data/genes_seed.txt --output results/both_results.txt
        """
    )
    
    parser.add_argument(
        '--algorithm',
        type=str,
        choices=['guild', 'diamond', 'both'],
        default='both',
        help='Algoritmo a ejecutar: guild, diamond, o both (default: both)'
    )
    
    parser.add_argument(
        '--network',
        type=str,
        required=True,
        help='Archivo de red de interacciones (formato: node1 node2 weight)'
    )
    
    parser.add_argument(
        '--seeds',
        type=str,
        required=True,
        help='Archivo con genes semilla (uno por línea)'
    )
    
    parser.add_argument(
        '--output',
        type=str,
        required=True,
        help='Archivo de salida para los resultados'
    )
    
    parser.add_argument(
        '--alpha',
        type=float,
        default=0.85,
        help='Parámetro alpha para GUILD RWR (default: 0.85)'
    )
    
    parser.add_argument(
        '--max-genes',
        type=int,
        default=200,
        help='Número máximo de genes para DIAMOnD (default: 200)'
    )
    
    parser.add_argument(
        '--network-type',
        type=str,
        choices=['auto', 'guild', 'diamond', 'string'],
        default='auto',
        help='Tipo de formato de red (default: auto)'
    )
    
    args = parser.parse_args()
    
    try:
        # Inicializar propagador
        propagator = NetworkPropagation(args.network, args.seeds)
        
        # Cargar red y semillas
        propagator.load_network(args.network_type)
        propagator.load_seeds()
        
        # Ejecutar algoritmos según lo solicitado
        if args.algorithm in ['guild', 'both']:
            logger.info("\n" + "="*60)
            logger.info("EJECUTANDO GUILD")
            logger.info("="*60)
            guild_results = propagator.guild_random_walk_restart(alpha=args.alpha)
            
            if args.algorithm == 'guild':
                propagator.save_results(guild_results, args.output, 'guild')
            else:
                guild_output = args.output.replace('.txt', '_guild.txt')
                propagator.save_results(guild_results, guild_output, 'guild')
        
        if args.algorithm in ['diamond', 'both']:
            logger.info("\n" + "="*60)
            logger.info("EJECUTANDO DIAMOnD")
            logger.info("="*60)
            diamond_results = propagator.diamond_algorithm(max_genes=args.max_genes)
            
            if args.algorithm == 'diamond':
                propagator.save_results(diamond_results, args.output, 'diamond')
            else:
                diamond_output = args.output.replace('.txt', '_diamond.txt')
                propagator.save_results(diamond_results, diamond_output, 'diamond')
        
        logger.info("\n" + "="*60)
        logger.info("PROCESO COMPLETADO EXITOSAMENTE")
        logger.info("="*60)
        
    except Exception as e:
        logger.error(f"Error en la ejecución: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()