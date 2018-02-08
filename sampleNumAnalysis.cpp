#include "sampleNumAnalysis.hpp"

using namespace llvm;

namespace sampleNumAnalysis {
    
    char SampleNumberAnalysis::ID = 0;
    static RegisterPass<SampleNumberAnalysis> X("sampleNumAnalysis", "Sample number analysis Pass", false, false);
    
    SampleNumberAnalysis::SampleNumberAnalysis() : FunctionPass(ID) {}
    
    string SampleNumberAnalysis::getBound(Value *bound) {
        
        if (isa<Instruction>(bound)) {
            
            Instruction *inst = cast<Instruction>(bound);
            
            switch (inst->getOpcode()) {
                case Instruction::Add:
                    return "(" + getBound(inst->getOperand(0)) + " + " + getBound(inst->getOperand(1)) + ")";
                    break;
                case Instruction::Sub:
                    return "(" + getBound(inst->getOperand(0)) + " - " + getBound(inst->getOperand(1)) + ")";;
                    break;
                case Instruction::Mul:
                    return "(" + getBound(inst->getOperand(0)) + " * " + getBound(inst->getOperand(1)) + ")";;
                    break;
                case Instruction::FDiv:
                case Instruction::SDiv:
                case Instruction::UDiv:
                    return "(" + getBound(inst->getOperand(0)) + " / " + getBound(inst->getOperand(1)) + ")";;
                    break;
                case Instruction::Load:
                    return inst->getOperand(0)->getName().str();
                    break;
                case Instruction::Alloca:
                    return inst->getName();
                    break;
                default:
                    break;
            }
        }
        else if (isa<ConstantInt>(bound)) {
            return to_string(dyn_cast<ConstantInt>(bound)->getValue().getSExtValue());
        }
        return "";
    }
    
    int SampleNumberAnalysis::getConstantStride(Value* expr) {
        
        if (isa<Instruction>(expr)) {
            Instruction *inst = dyn_cast<Instruction>(expr);
            
            switch (inst->getOpcode()) {
                case Instruction::Add:
                    if (isa<Constant>(inst->getOperand(1))) {
                        return dyn_cast<ConstantInt>(inst->getOperand(1))->getValue().getSExtValue();
                    }
                    break;
                default:
                    break;
            }
        }
        
        return 0;
    }

    int SampleNumberAnalysis::valueOfExpression(llvm::Value *expr, std::map<Value *, int> counter) {
        
        if (isa<Instruction>(expr)) {
            
            Instruction *inst = cast<Instruction>(expr);
            
            switch (inst->getOpcode()) {
                case Instruction::Add:
                    return  valueOfExpression(inst->getOperand(0), counter) + valueOfExpression(inst->getOperand(1), counter);
                    break;
                case Instruction::Sub:
                    return valueOfExpression(inst->getOperand(0), counter) - valueOfExpression(inst->getOperand(1), counter);
                    break;
                case Instruction::Mul:
                    return  valueOfExpression(inst->getOperand(0), counter) * valueOfExpression(inst->getOperand(1), counter);
                    break;
                case Instruction::FDiv:
                case Instruction::SDiv:
                case Instruction::UDiv:
                    return valueOfExpression(inst->getOperand(0), counter) / valueOfExpression(inst->getOperand(1), counter);
                    break;
                case Instruction::Load:
                    return counter[expr];
                    break;
                case Instruction::Alloca:
                    return counter[expr];
                    break;
                default:
                    break;
            }
        }
        else if (isa<ConstantInt>(expr)) {
            return dyn_cast<ConstantInt>(expr)->getValue().getSExtValue();
        }
        return 0;
        
    }
    
    void SampleNumberAnalysis::calculateSampleNum(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops) {
        
        if (LoopRefTree->L != NULL) {
            /* calculate sample number */
            
            vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops_New = loops;
            loops_New.push_back(LoopRefTree);
            
            bool allConstantBounds = true;
            for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator lit = loops_New.begin(), elit = loops_New.end(); lit != elit; ++lit) {
                for (unsigned long i = 0; i < (*lit)->LIS->IDV->size(); i++) {
                    std::string s_tmp = getBound((*(*lit)->LIS->LB)[i].first);
                    for (unsigned long si = 0; si < s_tmp.length(); si++) {
                        if (std::isdigit(s_tmp[si]) == 0) {
                            allConstantBounds = false;
                            break;
                        }
                    }
                    if (allConstantBounds == false) {
                        break;
                    }
                    
                    s_tmp = getBound((*(*lit)->LIS->LB)[i].second);
                    for (unsigned long si = 0; si < s_tmp.length(); si++) {
                        if (std::isdigit(s_tmp[si]) == 0) {
                            allConstantBounds = false;
                            break;
                        }
                    }
                    if (allConstantBounds == false) {
                        break;
                    }
                }
                if (allConstantBounds == false) {
                    break;
                }
            }
            
            
            if (allConstantBounds == true) {
                /* all constant bounds, calculate the number */

//                errs() << "all constants \n";
                
                uint64_t sampling_num = 1;
                for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator lit = loops_New.begin(), elit = loops_New.end(); lit != elit; ++lit) {
                    for (unsigned long i = 0; i < (*lit)->LIS->IDV->size(); i++) {
                        int stride = getConstantStride((*(*lit)->LIS->INC)[i]);
                        if (stride < 0) {
                            sampling_num *= ((uint64_t) stoi(getBound((*(*lit)->LIS->LB)[i].first)) - (uint64_t) stoi(getBound((*(*lit)->LIS->LB)[i].second))) / (uint64_t) (-stride);
                        } else {
                            sampling_num *= ((uint64_t) stoi(getBound((*(*lit)->LIS->LB)[i].second)) - (uint64_t) stoi(getBound((*(*lit)->LIS->LB)[i].first))) / (uint64_t) (stride);
                        }
                    } 
                }
                
                sampleNum[LoopRefTree] = sampling_num * RANDOM_REF_SAMPLING_RATE;
                
//                errs() << "here sample num " << sampling_num << "\n";
                
            } else {
                
                std::map<Value*, int> counter;
                
                vector<bool> lbConstant(loops_New.size());
                vector<bool> ubConstant(loops_New.size());
                
                /* check constant bounds */
                for (unsigned long i = 0; i < loops_New.size(); i++) {
                    lbConstant[i] = true;
                    ubConstant[i] = true;
                    
                    std::string s_tmp = getBound((*(loops_New[i])->LIS->LB)[0].first);
                    for (unsigned long si = 0; si < s_tmp.length(); si++) {
                        if (std::isdigit(s_tmp[si]) == 0) {
                            lbConstant[i] = false;
                            break;
                        }
                    }
                    
                    s_tmp = getBound((*(loops_New[i])->LIS->LB)[0].second);
                    for (unsigned long si = 0; si < s_tmp.length(); si++) {
                        if (std::isdigit(s_tmp[si]) == 0) {
                            ubConstant[i] = false;
                            break;
                        }
                    }
                }

                /* init counter */
                for (unsigned long i = 0; i < loops_New.size(); i++) {
                    if (i == 0) {
                        counter[(*loops_New[i]->LIS->IDV)[0]] = stoi(getBound((*loops_New[i]->LIS->LB)[0].first));
                    } else {
                        if (lbConstant[i] == true) {
                            counter[(*loops_New[i]->LIS->IDV)[0]] = stoi(getBound((*loops_New[i]->LIS->LB)[0].first));
                        } else {
                            counter[(*loops_New[i]->LIS->IDV)[0]] = valueOfExpression((*loops_New[i]->LIS->LB)[0].first, counter);
                        }
                    }
                    
//                    errs() << counter[(*loops_New[i]->LIS->IDV)[0]] << " ";
                    
                }
//                errs() << "\n";
                
                
                /* iterate to get sample number */
                
                uint64_t sampling_num = 0;
                int i_idx = loops_New.size() - 1;
                while(1) {
                    

//                    errs() << "counter: ";
//                    for (unsigned long i = 0; i < loops_New.size(); i++) {
//                        errs() << counter[(*loops_New[i]->LIS->IDV)[0]] << " ";
//                    }
//                    errs() << getBound((*loops_New[i_idx]->LIS->LB)[0].second);
//                    errs() << "\n";
                    
//                    (*loops_New[i_idx]->LIS->LB)[0].second->dump();
                    
                    int last_level_ub = valueOfExpression((*loops_New[i_idx]->LIS->LB)[0].second, counter);
                    
//                    errs() << "last level: " << last_level_ub << " - " << counter[(*loops_New[i_idx]->LIS->IDV)[0]] << "\n";
                    if (counter[(*loops_New[i_idx]->LIS->IDV)[0]] < last_level_ub) {
                        sampling_num += last_level_ub - counter[(*loops_New[i_idx]->LIS->IDV)[0]];
                    }
                    
                    while(i_idx - 1 >= 0) {
                        i_idx -= 1;
                        counter[(*loops_New[i_idx]->LIS->IDV)[0]] ++;
                        if (counter[(*loops_New[i_idx]->LIS->IDV)[0]] < valueOfExpression((*loops_New[i_idx]->LIS->LB)[0].second, counter)) {
                            break;
                        }
                    }
                    
                    if (i_idx == 0 && counter[(*loops_New[i_idx]->LIS->IDV)[0]] == valueOfExpression((*loops_New[i_idx]->LIS->LB)[0].second, counter)) {
                        break;
                    }
                    i_idx++;
                    
                    while(i_idx < loops_New.size()) {
                        counter[(*loops_New[i_idx]->LIS->IDV)[0]] = valueOfExpression((*loops_New[i_idx]->LIS->LB)[0].first, counter);
                        i_idx++;
                    }
                    i_idx = loops_New.size() - 1;
                }
                
//                for (unsigned long i = 0; i < loops_New.size(); i++) {
//                    errs() << counter[(*loops_New[i]->LIS->IDV)[0]] << " ";
//                }
//                errs() << "\n";
                
//                errs() << "sample number " << sampling_num << "\n";
                
                sampleNum[LoopRefTree] = sampling_num * RANDOM_REF_SAMPLING_RATE;

            }
        }
        
        if (LoopRefTree->next != NULL) {
            for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops_New = loops;
                if (LoopRefTree->L != NULL) {
                    loops_New.push_back(LoopRefTree);
                }
                calculateSampleNum(*it, loops_New);
            }
        }
        
        return;
    }
    
    void SampleNumberAnalysis::dumpSampleNum(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LTroot, std::string prefix) {
        
        if (LTroot->L != NULL) {
            if (sampleNum.find(LTroot) != sampleNum.end()) {
                errs() << prefix << "Sample number: " << sampleNum[LTroot] << "\n";
            } else {
                errs() << prefix << "Sample number: N/A \n";
            }
        }
        
        if (LTroot->next != NULL) {
            for (std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LTroot->next->begin(), eit = LTroot->next->end(); it != eit; ++it) {
                dumpSampleNum(*it, prefix + "--");
            }
        }
        
        return;
    }
    
    bool SampleNumberAnalysis::runOnFunction(Function &F) {
        
        errs() << " /* Start to analysis the number of samples\n";
        
        /* reading info from previous passes */
        loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree = getAnalysis<loopAnalysis::LoopIndvBoundAnalysis>().LoopRefTree;
        
        errs() << "calculating:\n";
        
        calculateSampleNum(LoopRefTree, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>());
        
        errs() << "Dump tree:\n";
        
        dumpSampleNum(LoopRefTree, "--");
        
        errs() << " End of sample analysis */\n";
        
        return false;
    }
    
    void SampleNumberAnalysis::getAnalysisUsage(AnalysisUsage &AU) const {
        AU.setPreservesAll();
        AU.addRequired<LoopInfoWrapperPass>();
        AU.addRequired<idxAnalysis::IndexAnalysis>();
        AU.addRequired<argAnalysis::ArgumentAnalysis>();
        AU.addRequired<gVarAnalysis::GlobalVariableAnalysis>();
        AU.addRequired<loopAnalysis::LoopIndvBoundAnalysis>();
        
        return;
    }
    
}
